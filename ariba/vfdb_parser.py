import re
import pyfastaq
import gzip
import pandas as pd
import subprocess
import os
from ariba import common

class Error (Exception): pass

name_regex = re.compile(r'^(?P<vfdb_id>V\S+) \((?P<name>.*?)\) (?P<description>.*\]) \[(?P<genus_etc>.*)\]$')
desc_regex = re.compile(r'^[^[]*\[.*\((?P<vf_id>V\S+)\) .*\]$')

class VfdbParser:
    def __init__(self, infile, outprefix):
        self.infile = infile
        self.outprefix = outprefix


    @classmethod
    def _fa_header_to_name_pieces(cls, fa_header):
        m = name_regex.search(fa_header)
        if m is None:
            return None
        else:
            return tuple([m.group(x) for x in ['vfdb_id', 'name', 'description', 'genus_etc']])


    @staticmethod
    def _fa_header_to_name_and_metadata(fa_header):
        name_data = VfdbParser._fa_header_to_name_pieces(fa_header)
        if name_data is None:
            return fa_header, '.'
        else:
            vfdb_id, name, description, genus_etc = name_data
            vf_id = VfdbParser._name_piece_description_to_vfid(description)
            return name.replace(' ', '_') + '.' + vfdb_id + '.' + genus_etc.replace(' ', '_'), vf_id

        
    @classmethod
    def _name_piece_description_to_vfid(cls, description):
        n = desc_regex.search(description)
        if n is None:
            return description
        else:
            return n.group('vf_id')

        
    @classmethod
    def _load_VFs_xls_metadata(cls, outprefix):
        filename = 'VFs.xls.gz'
        # get tmpdir as defined in ref_genes_getter.py
        outprefix = os.path.abspath(outprefix)
        tmpdir = outprefix + '.tmp.download'
        
        if not os.path.exists(tmpdir):
            try:
                os.mkdir(tmpdir)
            except:
                raise Error('Error mkdir ' + tmpdir)
            
        zipfile = os.path.join(tmpdir, filename)
        common.download_file('http://www.mgc.ac.cn/VFs/Down/' + filename, zipfile)
        
        retcode = subprocess.call('gunzip -t ' + zipfile, shell=True)
        if retcode != 0:
            raise Error("Error opening for reading gzipped file '" + filename + "'")

        with gzip.open(zipfile, 'r') as handle:
            vfdb_meta_df = pd.read_excel(handle, header = 1, usecols = ['VFID', 'VF_Name', 'VFcategory', 'Function', 'Mechanism'], keep_default_na = False)
        vfdb_meta_df['description'] = vfdb_meta_df['VF_Name']+' ('+vfdb_meta_df['VFcategory']+'). '+vfdb_meta_df['Function']+' '+vfdb_meta_df['Mechanism']
        vfid_dict = dict(zip(vfdb_meta_df['VFID'], vfdb_meta_df['description']))
        del(vfdb_meta_df)

        return vfid_dict

    
    @classmethod
    def _vfid_to_metadata(cls, vfid_dict, vf_id):
        if vf_id is not None:
            vf_metadata = vfid_dict.get(vf_id)
            if vf_metadata is None:
                return ''
            else:
                return vf_metadata
        else:
            return ''       
        

    def run(self):
        file_reader = pyfastaq.sequences.file_reader(self.infile)
        fa_out = pyfastaq.utils.open_file_write(self.outprefix + '.fa')
        tsv_out = pyfastaq.utils.open_file_write(self.outprefix + '.tsv')

        vfid_dict = VfdbParser._load_VFs_xls_metadata(self.outprefix)
        
        for seq in file_reader:
            original_id = seq.id
            seq.id, description = self._fa_header_to_name_and_metadata(seq.id)
            vf_metadata = VfdbParser._vfid_to_metadata(vfid_dict, description)
            if description == '.':
                seq.id = original_id.split()[0]
            print(seq.id, '1', '0', '.', '.', description + ': '+ vf_metadata + ' Original name: ' + original_id, sep='\t', file=tsv_out)
            print(seq, file=fa_out)

        pyfastaq.utils.close(fa_out)
        pyfastaq.utils.close(tsv_out)

