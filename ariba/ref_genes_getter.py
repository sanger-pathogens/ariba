class Error (Exception): pass

import os
import re
import shutil
import tarfile
import pyfastaq
import urllib.request
import time
import json
from ariba import common, card_record, vfdb_parser


allowed_ref_dbs = {
    'argannot',
    'card',
    'plasmidfinder',
    'resfinder',
    'srst2_argannot',
    'vfdb_core',
    'vfdb_full',
}

argannot_ref = '"ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes",\nGupta et al 2014, PMID: 24145532\n'


class RefGenesGetter:
    def __init__(self, ref_db, version=None, debug=False):
        if ref_db not in allowed_ref_dbs:
            raise Error('Error in RefGenesGetter. ref_db must be one of: ' + str(allowed_ref_dbs) + ', but I got "' + ref_db)
        self.ref_db=ref_db
        self.debug = debug
        self.genetic_code = 11
        self.max_download_attempts = 3
        self.sleep_time = 2
        self.version = version
        pyfastaq.sequences.genetic_code = self.genetic_code


    def _download_file(self, url, outfile):
        print('Downloading "', url, '" and saving as "', outfile, '" ...', end='', sep='', flush=True)
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


    def _get_card_versions(self, tmp_file):
        print('Getting available CARD versions')
        self._download_file('https://card.mcmaster.ca/download', tmp_file)
        p = re.compile(r'''href="(/download/.*?broad.*?v([0-9]+\.[0-9]+\.[0-9]+)\.tar\.gz)"''')
        versions = {}

        with open(tmp_file) as f:
            for line in f:
                got = p.findall(line)
                for match in got:
                    key = tuple([int(x) for x in match[1].split('.')])
                    versions[key] = 'https://card.mcmaster.ca' + match[0]

        if len(versions) == 0:
            raise Error('Error getting CARD versions. Cannot continue')

        print('Found versions:')

        for key, url in sorted(versions.items()):
            print('.'.join([str(x) for x in key]), url, sep='\t')

        os.unlink(tmp_file)
        return versions


    def _get_from_card(self, outprefix):
        outprefix = os.path.abspath(outprefix)
        tmpdir = outprefix + '.download'
        current_dir = os.getcwd()

        try:
            os.mkdir(tmpdir)
            os.chdir(tmpdir)
        except:
            raise Error('Error mkdir/chdir ' + tmpdir)

        versions = self._get_card_versions('download.html')
        if self.version is not None:
            key = tuple([int(x) for x in self.version.split('.')])
            if key not in versions:
                raise Error('Error! Did not find requested version ' + self.version)
        else:
            key = sorted(list(versions.keys()))[-1]
            self.version = '.'.join([str(x) for x in key])

        print('Getting version', self.version)
        card_tarball_url = versions[key]
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

        final_fasta = outprefix + '.fa'
        final_tsv = outprefix + '.tsv'
        log_file = outprefix + '.log'
        f_out_fa = pyfastaq.utils.open_file_write(final_fasta)
        f_out_tsv = pyfastaq.utils.open_file_write(final_tsv)
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

                print(fasta, file=f_out_fa)

                if gi == 'NA':
                    gene_or_not = '0'
                    variant_only = '0'
                elif len(data['snps']) == 0:
                    gene_or_not = '1'
                    variant_only = '0'
                else:
                    gene_or_not = '1'
                    variant_only = '1'

                print(fasta.id, gene_or_not, variant_only, '.', '.', data['ARO_name'], sep='\t', file=f_out_tsv)

                if len(data['snps']) == 0 and data['ARO_description'] != '':
                    print(fasta.id, gene_or_not, variant_only, '.', '.', data['ARO_description'], sep='\t', file=f_out_tsv)
                else:
                    for snp in data['snps']:
                        if data['ARO_description'] != '':
                            print(fasta.id, gene_or_not, variant_only, snp, '.', data['ARO_description'], sep='\t', file=f_out_tsv)


        pyfastaq.utils.close(f_out_fa)
        pyfastaq.utils.close(f_out_tsv)
        pyfastaq.utils.close(f_out_log)
        os.chdir(current_dir)
        if not self.debug:
            shutil.rmtree(tmpdir)

        print('Extracted data and written ARIBA input files\n')
        print('Finished. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
        print('You can use them with ARIBA like this:')
        print('ariba prepareref -f', final_fasta, '-m', final_tsv, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print('"The Comprehensive Antibiotic Resistance Database", McArthur et al 2013, PMID: 23650175')
        print('and in your methods say that version', self.version, 'of the database was used')


    def _get_from_resfinder(self, outprefix):
        outprefix = os.path.abspath(outprefix)
        final_fasta = outprefix + '.fa'
        final_tsv = outprefix + '.tsv'
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
        fout_fa = pyfastaq.utils.open_file_write(final_fasta)
        fout_tsv = pyfastaq.utils.open_file_write(final_tsv)
        used_names = {}

        for filename in os.listdir():
            if filename.endswith('.fsa'):
                print('   ', filename)
                file_reader = pyfastaq.sequences.file_reader(filename)
                for seq in file_reader:
                    try:
                        prefix, suffix = seq.id.split('_', maxsplit=1)
                        description = 'Original name: ' + seq.id
                        seq.id = prefix + '.' + suffix
                    except:
                        description = '.'

                    # names are not unique across the files
                    if seq.id in used_names:
                        used_names[seq.id] += 1
                        seq.id += '_' + str(used_names[seq.id])
                    else:
                        used_names[seq.id] = 1

                    print(seq, file=fout_fa)
                    print(seq.id, '1', '0', '.', '.', description, sep='\t', file=fout_tsv)

        pyfastaq.utils.close(fout_fa)
        pyfastaq.utils.close(fout_tsv)
        print('\nFinished combining files\n')
        os.chdir(current_dir)
        if not self.debug:
            shutil.rmtree(tmpdir)
        print('Finished. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
        print('You can use them with ARIBA like this:')
        print('ariba prepareref -f', final_fasta, '-m', final_tsv, 'output_directory\n')
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
        final_fasta = outprefix + '.fa'
        final_tsv = outprefix + '.tsv'

        seq_reader = pyfastaq.sequences.file_reader(genes_file)
        f_out_tsv = pyfastaq.utils.open_file_write(final_tsv)
        f_out_fa = pyfastaq.utils.open_file_write(final_fasta)

        for seq in seq_reader:
            original_id = seq.id
            seq.id = re.sub(r'\((.*)\)', r'\1.', seq.id.split()[0])
            print(seq, file=f_out_fa)
            print(seq.id, '1', '0', '.', '.', 'Original name: ' + original_id, sep='\t', file=f_out_tsv)


        pyfastaq.utils.close(f_out_tsv)
        pyfastaq.utils.close(f_out_fa)
        if not self.debug:
            shutil.rmtree(tmpdir)

        print('Finished. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
        print('You can use them with ARIBA like this:')
        print('ariba prepareref -f', final_fasta, '-m', final_tsv, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print(argannot_ref)


    def _get_from_plasmidfinder(self, outprefix):
        outprefix = os.path.abspath(outprefix)
        final_fasta = outprefix + '.fa'
        final_tsv = outprefix + '.tsv'
        tmpdir = outprefix + '.tmp.download'
        current_dir = os.getcwd()

        try:
            os.mkdir(tmpdir)
            os.chdir(tmpdir)
        except:
            raise Error('Error mkdir/chdir ' + tmpdir)

        zipfile = 'plasmidfinder.zip'
        cmd = 'curl -X POST --data "folder=plasmidfinder&filename=plasmidfinder.zip" -o ' + zipfile + ' https://cge.cbs.dtu.dk/cge/download_data.php'
        print('Downloading data with:', cmd, sep='\n')
        common.syscall(cmd)
        common.syscall('unzip ' + zipfile)

        print('Combining downloaded fasta files...')
        fout_fa = pyfastaq.utils.open_file_write(final_fasta)
        fout_tsv = pyfastaq.utils.open_file_write(final_tsv)
        name_count = {}

        for filename in os.listdir(tmpdir):
            if filename.endswith('.fsa'):
                print('   ', filename)
                file_reader = pyfastaq.sequences.file_reader(os.path.join(tmpdir, filename))
                for seq in file_reader:
                    original_id = seq.id
                    seq.id = seq.id.replace('_', '.', 1)
                    if seq.id in name_count:
                        name_count[seq.id] += 1
                        seq.id = seq.id + '.' + str(name_count[seq.id])
                    else:
                        name_count[seq.id] = 1

                    print(seq, file=fout_fa)
                    print(seq.id, '0', '0', '.', '.', 'Original name was ' + original_id, sep='\t', file=fout_tsv)

        pyfastaq.utils.close(fout_fa)
        pyfastaq.utils.close(fout_tsv)
        print('\nFinished combining files\n')
        os.chdir(current_dir)
        if not self.debug:
            shutil.rmtree(tmpdir)
        print('Finished. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
        print('You can use them with ARIBA like this:')
        print('ariba prepareref -f', final_fasta, '-m', final_tsv, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print('"PlasmidFinder and pMLST: in silico detection and typing of plasmids", Carattoli et al 2014, PMID: 24777092\n')


    def _get_from_srst2_argannot(self, outprefix):
        srst2_version = '0.2.0'
        srst2_url = 'https://github.com/katholt/srst2/raw/v' + srst2_version + '/data/ARGannot.r1.fasta'
        srst2_fa = outprefix + '.original.fa'
        command = 'wget -O ' + srst2_fa + ' ' + srst2_url
        common.syscall(command, verbose=True)

        final_fasta = outprefix + '.fa'
        final_tsv = outprefix + '.tsv'

        f_out_fa = pyfastaq.utils.open_file_write(final_fasta)
        f_out_meta = pyfastaq.utils.open_file_write(final_tsv)
        seq_reader = pyfastaq.sequences.file_reader(srst2_fa)

        for seq in seq_reader:
            original_id = seq.id
            name, extra = seq.id.split()
            cluster_id, cluster_name, allele_name, allele_id = name.split('__')
            seq.id = cluster_name + '.' + name
            print(seq, file=f_out_fa)
            print(seq.id, 1, 0, '.', '.', 'Original name: ' + original_id, sep='\t', file=f_out_meta)

        pyfastaq.utils.close(f_out_fa)
        pyfastaq.utils.close(f_out_meta)
        if not self.debug:
            os.unlink(srst2_fa)

        print('Finished downloading and converting data. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
        print('You can use them with ARIBA like this:')
        print('ariba prepareref -f', final_fasta, '-m', final_tsv, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print('"SRST2: Rapid genomic surveillance for public health and hospital microbiology labs",\nInouye et al 2014, Genome Medicine, PMID: 25422674\n')
        print(argannot_ref)
        print('and in your methods say that the ARG-ANNOT sequences were used from version', srst2_version, 'of SRST2.')


    def _get_from_vfdb_core(self, outprefix):
        self._get_from_vfdb_common(outprefix, 'VFDB_setA_nt.fas.gz','core')


    def _get_from_vfdb_full(self, outprefix):
         self._get_from_vfdb_common(outprefix, 'VFDB_setB_nt.fas.gz','full')


    def _get_from_vfdb_common(self, outprefix, filename, info_text):
        outprefix = os.path.abspath(outprefix)
        tmpdir = outprefix + '.tmp.download'

        try:
            os.mkdir(tmpdir)
        except:
            raise Error('Error mkdir ' + tmpdir)

        zipfile = os.path.join(tmpdir, filename)
        self._download_file('http://www.mgc.ac.cn/VFs/Down/' + filename, zipfile)
        print('Extracting files ... ', end='', flush=True)
        vparser = vfdb_parser.VfdbParser(zipfile, outprefix)
        vparser.run()
        if not self.debug:
            shutil.rmtree(tmpdir)
        print('done')
        final_fasta = outprefix + '.fa'
        final_tsv = outprefix + '.tsv'

        print('Extracted core DNA sequence dataset and metadata. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
        print('You can use them with ARIBA like this:')
        print('ariba prepareref -f', final_fasta, '-m', final_tsv, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print('"VFDB 2016: hierarchical and refined dataset for big data analysis-10 years on",\nChen LH et al 2016, Nucleic Acids Res. 44(Database issue):D694-D697. PMID: 26578559\n')


    def run(self, outprefix):
        exec('self._get_from_' + self.ref_db + '(outprefix)')

