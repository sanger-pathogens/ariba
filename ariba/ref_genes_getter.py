class Error (Exception): pass

import sys
import os
import shutil
import re
import requests
import pyfastaq
import urllib
import time
from bs4 import BeautifulSoup
from ariba import common


class RefGenesGetter:
    def __init__(self, ref_db, genetic_code=11):
        allowed_ref_dbs = {'card', 'argannot', 'resfinder'}
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


    def _get_souped_request(self, url):
        print('Getting url "', url, '" ...', sep='', end='')
        for i in range(self.max_download_attempts):
            time.sleep(self.sleep_time)
            r = requests.get(url)
            if r.status_code == 200:
                break
        else:
            raise Error('\nError requests.get with url: ' + url)

        print('done', flush=True)
        return BeautifulSoup(r.text, 'html.parser')


    def _get_card_gene_variant_info(self, gene, index_url):
        print('Getting variant info on CARD gene', gene, flush=True)
        soup = self._get_souped_request(index_url)

        # get link to Antibiotic Resistance page
        rows = soup.find_all('tr')
        gene_indexes = [i for i, j in enumerate(rows) if 'Antibiotic Resistance' in j.text]
        if len(gene_indexes) != 1:
            raise Error('Error getting one link to antibiotic resistance. Found ' + str(len(gene_indexes)) + ' links')

        row_index = gene_indexes[0] + 1
        assert row_index < len(rows)
        antibio_links = rows[row_index].find_all('a')
        print('Found', len(antibio_links), 'links to antibiotic resistance pages')

        variants = []

        for antibio_link_obj in antibio_links:
            antibio_link = antibio_link_obj['href']
            soup = self._get_souped_request(antibio_link)

            # get description
            ontology_def = soup.find(id='ontology-definition-field')
            if ontology_def is None:
                description = None
            else:
                # there are sometimew newline characters in the description.
                # replace all whitepace characters with a space
                description = re.sub('\s', ' ', ontology_def.text)

            # get variants
            bioinf_tables = [x for x in soup.find_all('table') if 'Bioinformatics' in x.text]

            if len(bioinf_tables) != 1:
                raise Error('Error getting Bioinformatics table from ' + antibio_link)

            bioinf_table = bioinf_tables[0]
            variant_elements = [x for x in bioinf_table.find_all('small') if 'Resistance Variant' in x.text]

            if len(variant_elements) < 1:
                print('WARNING:', gene, 'No variants found on page', antibio_link)
            else:
                new_variants = [x.text.split()[-1].split('<')[0] for x in variant_elements]
                for variant in new_variants:
                    print('New variant:', variant, description)
                    variants.append((variant, description))


        if len(variants):
            print('Total of', len(variants), 'variants found for gene', gene)
        else:
            print('WARNING:', gene, 'No valid variants found for gene')

        return variants


    def _get_card_variant_data(self, outfilehandle):
        soup = self._get_souped_request('http://arpcard.mcmaster.ca/?q=CARD/search/mqt.35950.mqt.806')
        table = soup.find(id='searchresultsTable')
        links = {x.text : x['href'] for x in table.find_all('a')}

        for gene, url in sorted(links.items()):
            print('\nGetting info for gene', gene, 'from', url)
            variants = self._get_card_gene_variant_info(gene, links[gene])
            if len(variants) == 0:
                print('WARNING: No valid variants found for gene', gene)
            for variant, description in variants:
                print(gene, 'p', variant, description, sep='\t', file=outfilehandle)


    @staticmethod
    def _card_parse_presence_absence(infile, outfile, metadata_fh):
        presence_absence_ids = set()

        file_reader = pyfastaq.sequences.file_reader(infile)
        f_out = pyfastaq.utils.open_file_write(outfile)

        for seq in file_reader:
            try:
                seq.id, description = seq.id.split(maxsplit=1)
            except:
                description = None

            presence_absence_ids.add(seq.id)
            print(seq, file=f_out)

            if description is not None:
                print(seq.id, '.', '.', description, sep='\t', file=metadata_fh)

        pyfastaq.utils.close(f_out)
        return presence_absence_ids


    @staticmethod
    def _card_parse_all_genes(infile, outfile, metadata_fh, presence_absence_ids):
        file_reader = pyfastaq.sequences.file_reader(infile)
        f_out = pyfastaq.utils.open_file_write(outfile)

        for seq in file_reader:
            try:
                seq.id, description = seq.id.split(maxsplit=1)
            except:
                description = None

            if seq.id in presence_absence_ids:
                continue

            print(seq, file=f_out)
            if description is not None:
                print(seq.id, '.', '.', description, sep='\t', file=metadata_fh)

        pyfastaq.utils.close(f_out)


    def _get_from_card(self, outprefix):
        tsv_fh = pyfastaq.utils.open_file_write(outprefix + '.metadata.tsv')
        self._get_card_variant_data(tsv_fh)

        print('Finished getting variant data. Getting FASTA files', flush=True)
        all_ref_genes_fa_gz = outprefix + '.tmp.downloaded.all_genes.fa.gz'
        presence_absence_fa_gz = outprefix + '.tmp.download.presence_absence.fa.gz'
        variants_only_fa = outprefix + '.variants_only.fa'
        presence_absence_fa = outprefix + '.presence_absence.fa'
        self._download_file('http://arpcard.mcmaster.ca/blast/db/nucleotide/AR-genes.fa.gz', all_ref_genes_fa_gz)
        self._download_file('http://arpcard.mcmaster.ca/blast/db/nucleotide/ARmeta-genes.fa.gz', presence_absence_fa_gz)

        print('Making presence_absence and variants_only fasta files, and udating metadata', flush=True)
        presence_absence_ids = self._card_parse_presence_absence(presence_absence_fa_gz, presence_absence_fa, tsv_fh)
        self._card_parse_all_genes(all_ref_genes_fa_gz, variants_only_fa, tsv_fh, presence_absence_ids)
        pyfastaq.utils.close(tsv_fh)

        print('Deleting temporary downloaded files')
        os.unlink(all_ref_genes_fa_gz)
        os.unlink(presence_absence_fa_gz)
        print('All done')


    def _get_from_resfinder(self, outprefix):
        outprefix = os.path.abspath(outprefix)
        final_fasta = outprefix + '.genes.fa'
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
                file_reader = pyfastaq.sequences.file_reader(os.path.join('database', filename))
                for seq in file_reader:
                    print(seq, file=f)

        pyfastaq.utils.close(f)

        print('\nCombined files. Final genes file is callled', final_fasta, end='\n\n')
        os.chdir(current_dir)
        shutil.rmtree(tmpdir)

        print('You can use it with ARIBA like this:')
        print('ariba run --presabs', os.path.relpath(final_fasta), 'reads_1.fq reads_2.fq output_directory\n')
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
        pyfastaq.tasks.to_fasta(genes_file, final_fasta)
        shutil.rmtree(tmpdir)

        print('Finished. Final genes file is called', final_fasta, end='\n\n')
        print('You can use it with ARIBA like this:')
        print('ariba run --presabs', os.path.relpath(final_fasta), 'reads_1.fq reads_2.fq output_directory\n')
        print('If you use this downloaded data, please cite:')
        print('"ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes",\nGupta et al 2014, PMID: 24145532\n')

    def run(self, outprefix):
        exec('self._get_from_' + self.ref_db + '(outprefix)')

