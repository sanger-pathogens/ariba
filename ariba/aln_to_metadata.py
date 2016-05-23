import os
import sys
import shutil
import pyfastaq

class Error (Exception): pass

class AlnToMetadata:
    def __init__(self,
      aln_file,
      vars_file,
    ):
        self.padded_seqs = AlnToMetadata._load_aln_file(aln_file)
        self.variants = AlnToMetadata._load_vars_file(vars_file)


    @classmethod
    def _load_aln_file(cls, aln_file):
        pass


    @classmethod
    def _load_vars_file(cls, vars_file):
        pass


    def run(self, outprefix):
        pass
