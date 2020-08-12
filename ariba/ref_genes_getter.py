class Error (Exception): pass

import os
import re
import tarfile
import pyfastaq
import json
import subprocess
import sys
from ariba import common, card_record, vfdb_parser, megares_data_finder, megares_zip_parser


allowed_ref_dbs = {
    'argannot',
    'card',
    'megares',
    'plasmidfinder',
    'resfinder',
    'srst2_argannot',
    'vfdb_core',
    'vfdb_full',
    'virulencefinder',
    'ncbi',#added by schultzm
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


    def _get_card_versions(self, tmp_file):
        print('Getting available CARD versions')
        common.download_file('https://card.mcmaster.ca/download', tmp_file, max_attempts=self.max_download_attempts, sleep_time=self.sleep_time, verbose=True)
        p = re.compile(r'''href="(/download/0/.*?broad.*?v([0-9]+\.[0-9]+\.[0-9]+)\.tar\.(gz|bz2))"''')
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
        card_tarball = 'card.tar.bz2'
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
            data['ARO_description'] = data['ARO_description'].encode('utf-8')
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
            common.rmtree(tmpdir)

        print('Extracted data and written ARIBA input files\n')
        print('Finished. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
        print('You can use them with ARIBA like this:')
        print('ariba prepareref -f', final_fasta, '-m', final_tsv, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print('"CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database", Alcock et al 2020, PMID: 31665441')
        print('and in your methods say that version', self.version, 'of the database was used')


    @classmethod
    def _get_genetic_epi_database_from_bitbucket(cls, db_name, outdir, git_commit=None):
        assert db_name in {'plasmidfinder', 'resfinder', 'virulencefinder'}
        cmd = 'git clone ' + 'https://bitbucket.org/genomicepidemiology/' + db_name + '_db.git ' + outdir
        common.syscall(cmd)

        if git_commit is not None:
            common.syscall('cd ' + outdir + ' && git checkout ' + git_commit)

        print('Using this git commit for ' + db_name + ' database:')
        subprocess.check_call('cd ' + outdir + ' && git log -n 1', shell=True)


    def _get_from_resfinder(self, outprefix):
        outprefix = os.path.abspath(outprefix)
        final_fasta = outprefix + '.fa'
        final_tsv = outprefix + '.tsv'
        tmpdir = outprefix + '.tmp.download'
        current_dir = os.getcwd()

        if self.version =='old':
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
        else:
            RefGenesGetter._get_genetic_epi_database_from_bitbucket('resfinder', tmpdir, git_commit=self.version)
            os.chdir(tmpdir)


        print('Combining downloaded fasta files...')
        fout_fa = pyfastaq.utils.open_file_write(final_fasta)
        fout_tsv = pyfastaq.utils.open_file_write(final_tsv)
        used_names = {}

        for filename in os.listdir():
            if filename.endswith('.fsa'):
                print('   ', filename)
                fix_file = os.path.join(tmpdir, filename + '.fix.fsa')
                RefGenesGetter._fix_virulencefinder_fasta_file(os.path.join(tmpdir, filename), fix_file)
                file_reader = pyfastaq.sequences.file_reader(fix_file)
                for seq in file_reader:
                    try:
                        prefix, suffix = seq.id.split('_', maxsplit=1)
                        description = 'Original name: ' + seq.id
                        seq.id = prefix + '.' + suffix
                    except:
                        description = '.'

                    # names are not unique across the files
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
            common.rmtree(tmpdir)
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

        downloadfile = 'arg-annot-database_doc'
        common.download_file('https://www.mediterranee-infection.com/wp-content/uploads/2019/09/ARG-ANNOT_NT_V6_July2019.txt', downloadfile, max_attempts=self.max_download_attempts, sleep_time=self.sleep_time, verbose=True)
        os.chdir(current_dir)
        print('Extracted files.')

        genes_file = os.path.join(tmpdir, downloadfile)
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
            common.rmtree(tmpdir)

        print('Finished. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
        print('You can use them with ARIBA like this:')
        print('ariba prepareref -f', final_fasta, '-m', final_tsv, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print(argannot_ref)


    def _get_from_megares(self, outprefix):
        data_finder = megares_data_finder.MegaresDataFinder(version=self.version)
        download_url = data_finder.run()
        zip_parser = megares_zip_parser.MegaresZipParser(download_url, outprefix)
        zip_parser.run()
        final_fasta = outprefix + '.fa'
        final_tsv = outprefix + '.tsv'
        print('Finished. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
        print('You can use them with ARIBA like this:')
        print('ariba prepareref -f', final_fasta, '-m', final_tsv, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print('"MEGARes: an antimicrobial database for high throughput sequencing", Lakin et al 2016, PMID: PMC5210519\n')


    def _get_from_plasmidfinder(self, outprefix):
        outprefix = os.path.abspath(outprefix)
        final_fasta = outprefix + '.fa'
        final_tsv = outprefix + '.tsv'
        tmpdir = outprefix + '.tmp.download'
        current_dir = os.getcwd()

        if self.version == 'old':
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
        else:
            RefGenesGetter._get_genetic_epi_database_from_bitbucket('plasmidfinder', tmpdir, git_commit=self.version)
            os.chdir(tmpdir)

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
            common.rmtree(tmpdir)
        print('Finished. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
        print('You can use them with ARIBA like this:')
        print('ariba prepareref -f', final_fasta, '-m', final_tsv, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print('"PlasmidFinder and pMLST: in silico detection and typing of plasmids", Carattoli et al 2014, PMID: 24777092\n')


    def _get_from_srst2_argannot(self, outprefix):
        if self.version is None:
            self.version = 'r2'
        if self.version not in {'r1', 'r2'}:
            raise Error('srst2_argannot version must be r1 or r2. Got this: ' + self.version)

        version_string = '.r1' if self.version == 'r1' else '_r2'
        srst2_url = 'https://raw.githubusercontent.com/katholt/srst2/master/data/ARGannot' + version_string + '.fasta'
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
        # Use to also output the version of SRST2 here, but the r2 version of their
        # fasta file was made after SRST2 release 0.2.0. At the time of writing this,
        # 0.2.0 is the latest release, ie r2 isn't in an SRST2 release.


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
        common.download_file('http://www.mgc.ac.cn/VFs/Down/' + filename, zipfile, max_attempts=self.max_download_attempts, sleep_time=self.sleep_time, verbose=True)
        print('Extracting files ... ', end='', flush=True)
        vparser = vfdb_parser.VfdbParser(zipfile, outprefix)
        vparser.run()
        if not self.debug:
            common.rmtree(tmpdir)
        print('done')
        final_fasta = outprefix + '.fa'
        final_tsv = outprefix + '.tsv'

        print('Extracted core DNA sequence dataset and metadata. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
        print('You can use them with ARIBA like this:')
        print('ariba prepareref -f', final_fasta, '-m', final_tsv, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print('"VFDB 2016: hierarchical and refined dataset for big data analysis-10 years on",\nChen LH et al 2016, Nucleic Acids Res. 44(Database issue):D694-D697. PMID: 26578559\n')


    @classmethod
    def _fix_virulencefinder_fasta_file(cls, infile, outfile):
        '''Some line breaks are missing in the FASTA files from
        virulence finder. Which means there are lines like this:
        AAGATCCAATAACTGAAGATGTTGAACAAACAATTCATAATATTTATGGTCAATATGCTATTTTCGTTGA
        AGGTGTTGCGCATTTACCTGGACATCTCTCTCCATTATTAAAAAAATTACTACTTAAATCTTTATAA>coa:1:BA000018.3
        ATGAAAAAGCAAATAATTTCGCTAGGCGCATTAGCAGTTGCATCTAGCTTATTTACATGGGATAACAAAG
        and therefore the sequences are messed up when we parse them. Also
        one has a > at the end, then the seq name on the next line.
        This function fixes the file by adding line breaks'''
        with open(infile) as f_in, open(outfile, 'w') as f_out:
            for line in f_in:
                if line.startswith('>') or '>' not in line:
                    print(line, end='', file=f_out)
                elif line.endswith('>\n'):
                    print('WARNING: found line with ">" at the end! Fixing. Line:' + line.rstrip() + ' in file ' + infile, file=sys.stderr)
                    print(line.rstrip('>\n'), file=f_out)
                    print('>', end='', file=f_out)
                else:
                    print('WARNING: found line with ">" not at the start! Fixing. Line:' + line.rstrip() + ' in file ' + infile, file=sys.stderr)
                    line1, line2 = line.split('>')
                    print(line1, file=f_out)
                    print('>', line2, sep='', end='', file=f_out)


    def _get_from_virulencefinder(self, outprefix):
        outprefix = os.path.abspath(outprefix)
        final_fasta = outprefix + '.fa'
        final_tsv = outprefix + '.tsv'
        tmpdir = outprefix + '.tmp.download'
        current_dir = os.getcwd()

        if self.version == 'old':
            try:
                os.mkdir(tmpdir)
                os.chdir(tmpdir)
            except:
                raise Error('Error mkdir/chdir ' + tmpdir)

            zipfile = 'virulencefinder.zip'
            cmd = 'curl -X POST --data "folder=virulencefinder&filename=virulencefinder.zip" -o ' + zipfile + ' https://cge.cbs.dtu.dk/cge/download_data.php'
            print('Downloading data with:', cmd, sep='\n')
            common.syscall(cmd)
            common.syscall('unzip ' + zipfile)
        else:
            RefGenesGetter._get_genetic_epi_database_from_bitbucket('virulencefinder', tmpdir, git_commit=self.version)
            os.chdir(tmpdir)

        print('Combining downloaded fasta files...')
        fout_fa = pyfastaq.utils.open_file_write(final_fasta)
        fout_tsv = pyfastaq.utils.open_file_write(final_tsv)
        name_count = {}

        for filename in os.listdir(tmpdir):
            if filename.endswith('.fsa'):
                print('   ', filename)
                fix_file = os.path.join(tmpdir, filename + '.fix.fsa')
                RefGenesGetter._fix_virulencefinder_fasta_file(os.path.join(tmpdir, filename), fix_file)
                file_reader = pyfastaq.sequences.file_reader(fix_file)
                for seq in file_reader:
                    original_id = seq.id
                    seq.id = seq.id.replace('_', '.', 1)
                    seq.id = seq.id.replace(' ', '_')
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
            common.rmtree(tmpdir)
        print('Finished. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
        print('You can use them with ARIBA like this:')
        print('ariba prepareref -f', final_fasta, '-m', final_tsv, 'output_directory\n')
        print('If you use this downloaded data, please cite:')
        print('"Real-time whole-genome sequencing for routine typing, surveillance, and outbreak detection of verotoxigenic Escherichia coli", Joensen al 2014, PMID: 24574290\n')



    def _get_from_ncbi(self, outprefix, test=None): ## author github:schultzm
        """
        Download the NCBI-curated Bacterial Antimicrobial Resistance Reference Gene Database.
        Uses BioPython to do the data collection and extraction.
        Author: schultzm (github) May 31, 2019.

        >>> from Bio import Entrez
        >>> import getpass
        >>> import socket
        >>> BIOPROJECT = "PRJNA313047"
        >>> RETMAX = 100
        >>> import getpass
        >>> import socket
        >>> Entrez.email = getpass.getuser()+'@'+socket.getfqdn()
        >>> search_results = Entrez.read(Entrez.esearch(db="nucleotide",
        ...                                             term=BIOPROJECT,
        ...                                             retmax=RETMAX,
        ...                                             usehistory="y",
        ...                                             idtype="acc"))
        >>> webenv = search_results["WebEnv"]
        >>> query_key = search_results["QueryKey"]
        >>> target_accn = "    NG_061627.1"
        >>> records = Entrez.efetch(db="nucleotide",
        ...                         rettype="gbwithparts",
        ...                         retmode="text",
        ...                         retstart=0,
        ...                         retmax=RETMAX,
        ...                         webenv=webenv,
        ...                         query_key=query_key,
        ...                         idtype="acc")
        >>> from Bio.Alphabet import generic_dna
        >>> from Bio import SeqIO
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> for gb_record in SeqIO.parse(records, "genbank"):
        ...     if gb_record.id == 'NG_061627.1':
        ...         gb_record
        SeqRecord(seq=Seq('TAATCCTTGGAAACCTTAGAAATTGATGGAGGATCTTAACAAGATCCTGACATA...GGC', IUPACAmbiguousDNA()), id='NG_061627.1', name='NG_061627', description='Klebsiella pneumoniae SCKLB88 mcr-8 gene for phosphoethanolamine--lipid A transferase MCR-8.2, complete CDS', dbxrefs=['BioProject:PRJNA313047'])

        """

        outprefix = os.path.abspath(outprefix)
        final_fasta = outprefix + '.fa'
        final_tsv = outprefix + '.tsv'

        # Download the database as genbank using Bio.Entrez
        from Bio import Entrez
        import getpass
        import socket
        import sys
        BIOPROJECT = "PRJNA313047" ## https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA313047
        RETMAX=100000000
        Entrez.email = getpass.getuser()+'@'+socket.getfqdn()
        # See section 9.15  Using the history and WebEnv in
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:entrez-webenv
        search_results = Entrez.read(Entrez.esearch(db="nucleotide",
                                                    term=BIOPROJECT,
                                                    retmax=RETMAX,
                                                    usehistory="y",
                                                    idtype="acc"))
        acc_list = search_results["IdList"]
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]
        if test:
            return acc_list
            #up to here
        if len(acc_list) > 0:
            print(f"E-fetching {len(acc_list)} genbank records from BioProject {BIOPROJECT} and writing to.  This may take a while.", file=sys.stderr)
            records = Entrez.efetch(db="nucleotide",
                                               rettype="gbwithparts", retmode="text",
                                               retstart=0, retmax=RETMAX,
                                               webenv=webenv, query_key=query_key,
                                               idtype="acc")
            #pull out the records as fasta from the genbank
            from Bio.Alphabet import generic_dna
            from Bio import SeqIO
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord

            print(f"Parsing genbank records.")
            with open(final_fasta, "w") as f_out_fa, \
                 open(final_tsv, "w") as f_out_tsv:
                for idx, gb_record in enumerate(SeqIO.parse(records, "genbank")):
                    print(f"'{gb_record.id}'")
                    n=0
                    record_new=[]
                    for index, feature in enumerate(gb_record.features):
                        if feature.type == 'CDS':
                            n+=1
                            gb_feature = gb_record.features[index]
                            id = None
                            try:
                                id = gb_feature.qualifiers["allele"]
                            except:
                                try:
                                    try:
                                        id = gb_feature.qualifiers["gene"]
                                    except:
                                        id = gb_feature.qualifiers["locus_tag"]
                                except KeyError:
                                    print(f"gb_feature.qualifer not found", file=sys.stderr)
                            accession = gb_record.id
                            seq_out = Seq(str(gb_feature.extract(gb_record.seq)), generic_dna)
                            record_new.append(SeqRecord(seq_out,
                                         id=f"{id[0]}.{accession}",
                                         description=""))
                    if len(record_new) == 1:
                        print(f"Processing record {idx+1} of {len(acc_list)} (accession {accession})", file=sys.stderr)
                        f_out_fa.write(f"{record_new[0].format('fasta').rstrip()}\n")
                        f_out_tsv.write(f"{id[0]}.{accession}\t1\t0\t.\t.\t{gb_feature.qualifiers['product'][0]}\n")
                    if idx == len(acc_list)-1:
                        print('Finished. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
                        print('You can use them with ARIBA like this:')
                        print('ariba prepareref -f', final_fasta, '-m', final_tsv, 'output_directory\n')

        else:
            print(f"Nothing to do. Exiting.")
    def run(self, outprefix):
        exec('self._get_from_' + self.ref_db + '(outprefix)')
