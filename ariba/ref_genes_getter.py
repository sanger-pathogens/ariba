class Error (Exception): pass

import os
import sys
import shutil
import tarfile
import pyfastaq
import urllib.request
import time
import json
from ariba import common, card_record, vfdb_parser


class RefGenesGetter:
    def __init__(self, ref_db, genetic_code=11):
        allowed_ref_dbs = {'card', 'argannot', 'resfinder','vfdb'}
        if ref_db not in allowed_ref_dbs:
            raise Error('Error in RefGenesGetter. ref_db must be one of: ' + str(allowed_ref_dbs) + ', but I got "' + ref_db)
        self.ref_db=ref_db
        self.genetic_code = genetic_code
        self.max_download_attempts = 3
        self.sleep_time = 2
        pyfastaq.sequences.genetic_code = self.genetic_code


    def _download_file(self, url, outfile):
        print('Downloading "', url, '" and saving as "', outfile, '" ...', end='', sep='')
        for i in range(self.max_download_attempts):
            time.sleep(self.sleep_time)
            try:
                urllib.request.urlretrieve(url, filename=outfile)
            except:
                continue
            break
        else:
            raise Error('Error downloading: ' + url)
        print(' done', flush=True)



    def _get_from_card(self, outprefix):
        outprefix = os.path.abspath(outprefix)
        tmpdir = outprefix + '.download'
        current_dir = os.getcwd()

        try:
            os.mkdir(tmpdir)
            os.chdir(tmpdir)
        except:
            raise Error('Error mkdir/chdir ' + tmpdir)

        card_version = '1.0.6'
        card_tarball_url = 'https://card.mcmaster.ca/download/0/broadsteet-v' + card_version + '.tar.gz'
        card_tarball = 'card.tar.gz'
        print('Working in temporary directory', tmpdir)
        print('Downloading data from card:', card_tarball_url, flush=True)
        common.syscall('wget -O ' + card_tarball + ' ' + card_tarball_url, verbose=True)
        print('...finished downloading', flush=True)
        if not tarfile.is_tarfile(card_tarball):
            raise Error('File ' + card_tarball + ' downloaded from ' + card_tarball_url + ' does not look like a valid tar archive. Cannot continue')

        json_file = './card.json'
        with tarfile.open(card_tarball, 'r') as tfile:
            tfile.extract(json_file)

        print('Extracted json data file ', json_file,'. Reading its contents...', sep='')

        variant_metadata_tsv = outprefix + '.metadata.tsv'
        presence_absence_fa = outprefix + '.presence_absence.fa'
        variants_only_fa = outprefix + '.variants_only.fa'
        noncoding_fa = outprefix + '.noncoding.fa'
        log_file = outprefix + '.log'
        f_out_tsv = pyfastaq.utils.open_file_write(variant_metadata_tsv)
        f_out_presabs = pyfastaq.utils.open_file_write(presence_absence_fa)
        f_out_var_only = pyfastaq.utils.open_file_write(variants_only_fa)
        f_out_noncoding = pyfastaq.utils.open_file_write(noncoding_fa)
        f_out_log = pyfastaq.utils.open_file_write(log_file)

        with open(json_file) as f:
            json_data = json.load(f)

        json_data = {int(x): json_data[x] for x in json_data if not x.startswith('_')}
        print('Found', len(json_data), 'records in the json file. Analysing...', flush=True)

        for gene_key, gene_dict in sorted(json_data.items()):
            crecord = card_record.CardRecord(gene_dict)
            data = crecord.get_data()
            fasta_name_prefix = '.'.join([
                card_record.CardRecord._ARO_name_to_fasta_name(data['ARO_name']),
                data['ARO_accession'],
            ])

            for card_key, gi, genbank_id, start, end, dna_seq, protein_seq in data['dna_seqs_and_ids']:
                if dna_seq == '':
                    print('Empty dna sequence', gene_key, data['ARO_id'], data['ARO_accession'], sep='\t', file=f_out_log)
                    continue

                fasta_id = '.'.join([
                    fasta_name_prefix,
                    genbank_id,
                    start + '-' + end,
                    card_key
                ])
                fasta = pyfastaq.sequences.Fasta(fasta_id, dna_seq)
                variant_type = 'p'

                if gi != 'NA':
                    gene_tuple = fasta.make_into_gene()
                    if gene_tuple is None:
                        print('Could not make gene from sequence', fasta.id, sep='\t', file=f_out_log)
                        continue
                    else:
                        translated =  gene_tuple[0].translate()
                        if gene_tuple[0][:3] in pyfastaq.genetic_codes.starts[self.genetic_code]:
                            translated.seq = 'M' + translated.seq[1:]

                        if translated.seq[:-1] != protein_seq:
                            print('Translation of inferred gene dna sequence does not match protein sequence', fasta.id, sep='\t', file=f_out_log)
                            continue

                if gi == 'NA':
                    fasta_filehandle = f_out_noncoding
                    variant_type = 'n'
                elif len(data['snps']) == 0:
                    fasta_filehandle = f_out_presabs
                else:
                    fasta_filehandle = f_out_var_only

                print(fasta.id, '.', '.', '.', data['ARO_name'], sep='\t', file=f_out_tsv)

                if len(data['snps']) == 0:
                    print(fasta, file=fasta_filehandle)
                    print(fasta.id, '.', '.', '.', data['ARO_description'], sep='\t', file=f_out_tsv)
                else:
                    print(fasta, file=fasta_filehandle)
                    for snp in data['snps']:
                        print(fasta.id, variant_type, snp, '.', data['ARO_description'], sep='\t', file=f_out_tsv)


        pyfastaq.utils.close(f_out_tsv)
        pyfastaq.utils.close(f_out_presabs)
        pyfastaq.utils.close(f_out_var_only)
        pyfastaq.utils.close(f_out_noncoding)
        pyfastaq.utils.close(f_out_log)
        os.chdir(current_dir)
        print('Extracted data and written ARIBA input files\n')
        print('Final genes files and metadata file:')
        print('   ', presence_absence_fa)
        print('   ', variants_only_fa)
        print('   ', variant_metadata_tsv)

        print('\nYou can use those files with ARIBA like this:')
        print('ariba prepareref --ref_prefix', outprefix, 'output_directory\n')

        print('If you use this downloaded data, please cite:')
        print('"The Comprehensive Antibiotic Resistance Database", McArthur et al 2013, PMID: 23650175')
        print('and in your methods say that version', card_version, 'of the database was used')


    def _get_from_resfinder(self, outprefix):
        outprefix = os.path.abspath(outprefix)
        final_fasta = outprefix + '.presence_absence.fa'
        tmpdir = outprefix + '.tmp.download'
        current_dir = os.getcwd()

        try:
            os.mkdir(tmpdir)
            os.chdir(tmpdir)
        except:
            raise Error('Error mkdir/chdir ' + tmpdir)

        zipfile = 'resfinder.zip'
        cmd = 'curl -X POST --data "folder=resfinder&filename=resfinder.zip" -o ' + zipfile + ' https://cge.cbs.dtu.dk/cge/download_data.php'
        print('Downloading data with:', cmd, sep='\n')
        common.syscall(cmd)
        common.syscall('unzip ' + zipfile)

        print('Combining downloaded fasta files...')
        f = pyfastaq.utils.open_file_write(final_fasta)

        for filename in os.listdir('database'):
            if filename.endswith('.fsa'):
                print('   ', filename)
                prefix = filename.split('.')[0]
                file_reader = pyfastaq.sequences.file_reader(os.path.join('database', filename))
                for seq in file_reader:
                    seq.id = prefix + '.' + seq.id
                    print(seq, file=f)

        pyfastaq.utils.close(f)

        print('\nCombined files. Final genes file is called', final_fasta, end='\n\n')
        os.chdir(current_dir)
        shutil.rmtree(tmpdir)

        print('You can use it with ARIBA like this:')
        print('ariba prepareref --ref_prefix', outprefix, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print('"Identification of acquired antimicrobial resistance genes", Zankari et al 2012, PMID: 22782487\n')


    def _get_from_argannot(self, outprefix):
        outprefix = os.path.abspath(outprefix)
        tmpdir = outprefix + '.tmp.download'
        current_dir = os.getcwd()

        try:
            os.mkdir(tmpdir)
            os.chdir(tmpdir)
        except:
            raise Error('Error mkdir/chdir ' + tmpdir)

        zipfile = 'arg-annot-database_doc.zip'
        self._download_file('http://www.mediterranee-infection.com/arkotheque/client/ihumed/_depot_arko/articles/304/arg-annot-database_doc.zip', zipfile)
        common.syscall('unzip ' + zipfile)
        os.chdir(current_dir)
        print('Extracted files.')

        genes_file = os.path.join(tmpdir, 'Database Nt Sequences File.txt')
        final_fasta = outprefix + '.presence_absence.fa'

        seq_reader = pyfastaq.sequences.file_reader(genes_file)
        ids = {}
        for seq in seq_reader:
            ids[seq.id] = ids.get(seq.id, 0) + 1

        for name, count in sorted(ids.items()):
            if count > 1:
                print('Warning! Sequence name', name, 'found', count, 'times in download. Keeping longest sequence', file=sys.stderr)

        pyfastaq.tasks.to_unique_by_id(genes_file, final_fasta)
        shutil.rmtree(tmpdir)

        print('Finished. Final genes file is called', final_fasta, end='\n\n')
        print('You can use it with ARIBA like this:')
        print('ariba prepareref --ref_prefix', outprefix, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print('"ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes",\nGupta et al 2014, PMID: 24145532\n')

    def _get_from_vfdb(self, outprefix):
        outprefix = os.path.abspath(outprefix)
        tmpdir = outprefix + '.tmp.download'
        current_dir = os.getcwd()

        try:
            os.mkdir(tmpdir)
        except:
            raise Error('Error mkdir ' + tmpdir)

        zipfile = os.path.join(tmpdir, 'VFDB_setA_nt.fas.gz')
        self._download_file('http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz', zipfile)
        vparser = vfdb_parser.VfdbParser(zipfile, outprefix)
        vparser.run()
        shutil.rmtree(tmpdir)
        print('Extracted files.')
        final_fasta = outprefix + '.presence_absence.fa'
        final_tsv = outprefix + '.metadata.tsv'

        print('Extracted core DNA sequence dataset and metadata. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
        print('You can use it with ARIBA like this:')
        print('ariba prepareref --ref_prefix', outprefix, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print('"VFDB 2016: hierarchical and refined dataset for big data analysis-10 years on",\nChen LH et al 2016, Nucleic Acids Res. 44(Database issue):D694-D697. PMID: 26578559\n')

    def run(self, outprefix):
        exec('self._get_from_' + self.ref_db + '(outprefix)')

