class Error (Exception): pass

import sys
import os
import shutil
import requests
import pyfastaq
import urllib
import time
from bs4 import BeautifulSoup
from ariba import refcheck, common


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

        variants = set()

        for antibio_link_obj in antibio_links:
            antibio_link = antibio_link_obj['href']
            soup = self._get_souped_request(antibio_link)

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
                variants.update(new_variants)


        if len(variants):
            print('Total of', len(variants), 'variants found for gene', gene)
        else:
            print('WARNING:', gene, 'No valid variants found for gene')

        return variants


    def _write_card_variants_tsv(self, variants, outfile):
        f = pyfastaq.utils.open_file_write(outfile)
        for gene in sorted(variants):
            print(gene,  ','.join(variants[gene]), sep='\t', file=f)
        pyfastaq.utils.close(f)


    def _get_from_card(self, outprefix):
        all_ref_genes_fa_gz = outprefix + '.all_ref.genes.00.fa.gz'
        self._download_file('http://arpcard.mcmaster.ca/blast/db/nucleotide/AR-genes.fa.gz', all_ref_genes_fa_gz)
        soup = self._get_souped_request('http://arpcard.mcmaster.ca/?q=CARD/search/mqt.35950.mqt.806')
        table=soup.find(id='searchresultsTable')
        links = {x.text : x['href'] for x in table.find_all('a')}
        variants = {}

        for gene, url in sorted(links.items()):
            print('\nGetting info for gene', gene, 'from', url, flush=True)
            variants_strings = self._get_card_gene_variant_info(gene, links[gene])
            if len(variants_strings) == 0:
                print('WARNING: No valid variants found for gene', gene)
            variants[gene] = variants_strings

        self._write_card_variants_tsv(variants, outprefix + '.variants.tsv')

        all_ref_genes_fa = outprefix + '.original_downloaded_genes.fa'
        pyfastaq.utils.syscall('gunzip -c ' + all_ref_genes_fa_gz + ' > ' + all_ref_genes_fa)
        print('\nProcessing reference gene fasta file', all_ref_genes_fa, flush=True)
        all_ref_genes_fa_fixnames = outprefix + '.tmp.all_ref.genes.01.fixnames.fa'
        pyfastaq.tasks.to_fasta(all_ref_genes_fa, all_ref_genes_fa_fixnames, strip_after_first_whitespace=True, check_unique=True)
        refcheck_prefix = outprefix + '.tmp.all_ref.genes.02.refcheck'
        refchecker = refcheck.Checker(all_ref_genes_fa_fixnames, outprefix=refcheck_prefix, genetic_code=self.genetic_code)
        refchecker.run()
        os.rename(refcheck_prefix + '.fa', outprefix + '.genes.fa')
        print('Finished processing reference gene fasta file')


    def _get_from_resfinder(self, outprefix):
        outprefix = os.path.abspath(outprefix)
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

        print('Combining downloaded fasta files')
        all_genes_fasta = outprefix + '.original_downloaded_genes.fa'
        f = pyfastaq.utils.open_file_write(all_genes_fasta)

        for filename in os.listdir():
            if filename.endswith('.fsa'):
                print('   ', filename)
                file_reader = pyfastaq.sequences.file_reader(filename)
                for seq in file_reader:
                    print(seq, file=f)

        pyfastaq.utils.close(f)

        print('Combined files. Running refcheck on genes file')
        refcheck_prefix = outprefix + '.tmp.refcheck'
        refchecker = refcheck.Checker(all_genes_fasta, outprefix=refcheck_prefix, genetic_code=self.genetic_code)
        refchecker.run()
        os.chdir(current_dir)
        shutil.rmtree(tmpdir)
        final_fasta = outprefix + '.genes.fa'
        os.rename(refcheck_prefix + '.fa', final_fasta)
        print('Finished. Final genes file is called', final_fasta)


    def _get_from_argannot(self, outprefix):
        refcheck_prefix = os.path.abspath(outprefix) + '.tmp.refcheck'
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
        genes_file = 'Database Nt Sequences File.txt'

        print('Extracted files. Running refcheck on genes file')
        refchecker = refcheck.Checker(genes_file, outprefix=refcheck_prefix, genetic_code=self.genetic_code)
        refchecker.run()
        os.chdir(current_dir)
        os.rename(os.path.join(tmpdir, genes_file), outprefix + '.original_downloaded_genes.fa')
        shutil.rmtree(tmpdir)
        final_fasta = outprefix + '.genes.fa'
        os.rename(refcheck_prefix + '.fa', final_fasta)
        print('Finished. Final genes file is called', final_fasta)


    def run(self, outprefix):
        exec('self._get_from_' + self.ref_db + '(outprefix)')

