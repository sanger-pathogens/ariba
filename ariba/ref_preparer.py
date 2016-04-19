import os
import sys
import pyfastaq
from ariba import cdhit, reference_data

class Error (Exception): pass


class RefPreparer:
    def __init__(self,
        ref_prefix=None,
        presabs=None,
        varonly=None,
        noncoding=None,
        metadata=None,
        min_gene_length=6,
        max_gene_length=10000,
        genetic_code=11,
        verbose=False,
    ):
        self.ref_prefix = ref_prefix
        self.presabs = presabs
        self.varonly = varonly
        self.noncoding = noncoding
        self.metadata = metadata
        self.min_gene_length = min_gene_length
        self.max_gene_length = max_gene_length
        self.genetic_code = genetic_code
        self.verbose = verbose


    @staticmethod
    def _get_ref_files(ref_prefix, presabs, varonly, noncoding, metadata, verbose=False):
        if {None} == {ref_prefix, presabs, varonly, noncoding}:
            raise Error('Error in RefPreparer._get_ref_files. All input files and ref_prefix were None. Cannot continue')

        filenames = {
            'presabs': presabs,
            'varonly': varonly,
            'noncoding': noncoding,
            'metadata': metadata,
        }

        file_suffixes = {
            'presabs': 'presence_absence.fa',
            'varonly': 'variants_only.fa',
            'noncoding': 'noncoding.fa',
            'metadata': 'metadata.tsv',
        }

        if verbose:
            print('Looking for input files ...')

        for key in file_suffixes:
            if ref_prefix is not None:
                filename = os.path.abspath(ref_prefix + '.' + file_suffixes[key])

                if os.path.exists(filename):
                    if verbose:
                        print('Found: ', filename, '.\n    ...treating it as if this was used: --', key, ' ', filename, sep='')
                    filenames[key] = filename
                else:
                    if verbose:
                        print('Not found:', filename)
                    filenames[key] = None
            elif filenames[key] is not None:
                if os.path.exists(filenames[key]):
                    filenames[key] = os.path.abspath(filenames[key])
                    if verbose:
                        print('Found: ', filename, 'from option --', key, sep='')
                else:
                    raise Error('File not found! Cannot continue. Looked for: ' + filenames[key])

        if {None} == {filenames['presabs'], filenames['varonly'], filenames['noncoding']}:
            raise Error('Error in RefPreparer._get_ref_files. No FASTA files given! Cannot continue')

        return filenames


    def run(self, outdir):
        if os.path.exists(outdir):
            raise Error('Error! Output directory ' + outdir + ' already exists. Cannot continue')

        try:
            os.mkdir(outdir)
        except:
            raise Error('Error making output directory ' + outdir + '. Cannot continue')


        filenames = self._get_ref_files(self.ref_prefix, self.presabs, self.varonly, self.metadata, self.verbose)

        refdata = reference_data.ReferenceData(
            presence_absence_fa=filenames['presabs'],
            variants_only_fa=filenames['varonly'],
            non_coding_fa=filenames['noncoding'],
            metadata_tsv=filenames['metadata'],
            min_gene_length=self.min_gene_length,
            max_gene_length=self.max_gene_length,
            genetic_code=self.genetic_code,
        )

        if self.verbose:
            print('{:_^79}'.format(' Checking reference data '), flush=True)
        refdata.sanity_check(os.path.join(outdir, 'refcheck'))

