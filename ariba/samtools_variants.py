import os
import sys
import pysam
import pyfastaq
from ariba import common

class Error (Exception): pass


class SamtoolsVariants:
    def __init__(self,
      ref_fa,
      bam,
      outprefix,
      log_fh=sys.stdout,
      samtools_exe='samtools',
      bcftools_exe='bcftools',
      bcf_min_dp=10,
      bcf_min_dv=5,
      bcf_min_dv_over_dp=0.3,
      bcf_min_qual=20,
    ):
        self.ref_fa = os.path.abspath(ref_fa)
        self.bam = os.path.abspath(bam)
        self.outprefix = os.path.abspath(outprefix)
        self.log_fh = log_fh
        self.samtools_exe = samtools_exe
        self.bcftools_exe = bcftools_exe
        self.bcf_min_dp = bcf_min_dp
        self.bcf_min_dv = bcf_min_dv
        self.bcf_min_dv_over_dp = bcf_min_dv_over_dp
        self.bcf_min_qual = bcf_min_qual

        self.vcf_file = self.outprefix + '.vcf'
        self.read_depths_file = self.outprefix + '.read_depths.gz'


    def _make_vcf_and_read_depths_files(self):
        tmp_vcf = self.vcf_file + '.tmp'
        cmd = ' '.join([
            self.samtools_exe, 'mpileup',
            '-t INFO/AD',
            '-A',
            '-f', self.ref_fa,
            '-u',
            '-v',
            self.bam,
            '>',
            tmp_vcf
        ])

        common.syscall(cmd, verbose=True, verbose_filehandle=self.log_fh)

        cmd = ' '.join([
            self.bcftools_exe, 'call -m',
            tmp_vcf,
            '|',
            self.bcftools_exe, 'query',
            r'''-f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AD]\n' ''',
            '>',
            self.read_depths_file + '.tmp'
        ])

        common.syscall(cmd, verbose=True, verbose_filehandle=self.log_fh)
        pysam.tabix_compress(self.read_depths_file + '.tmp', self.read_depths_file)
        pysam.tabix_index(self.read_depths_file, seq_col=0, start_col=1, end_col=1)
        os.unlink(self.read_depths_file + '.tmp')

        cmd = ' '.join([
            self.bcftools_exe, 'call -m -v',
            tmp_vcf,
            '|',
            self.bcftools_exe, 'filter',
            '-i', '"SUM(AD)>=5 & MIN(AD)/DP>=0.1"',
            '-o', self.vcf_file
        ])

        common.syscall(cmd, verbose=True, verbose_filehandle=self.log_fh)
        os.unlink(tmp_vcf)


    @classmethod
    def _get_read_depths(cls, read_depths_file, sequence_name, position):
        '''Returns total read depth and depth of reads supporting alternative (if present)'''
        assert os.path.exists(read_depths_file)
        assert os.path.exists(read_depths_file + '.tbi')
        tbx = pysam.TabixFile(read_depths_file)
        try:
            rows = [x for x in tbx.fetch(sequence_name, position, position + 1)]
        except:
            return None

        if len(rows) > 1: # which happens with indels, mutiple lines for same base of reference
            test_rows = [x for x in rows if x.rstrip().split()[3] != '.']
            if len(test_rows) != 1:
                rows = [rows[-1]]
            else:
                rows = test_rows

        if len(rows) == 1:
            r, p, ref_base, alt_base, ref_counts, alt_counts = rows[0].rstrip().split()
            return ref_base, alt_base, int(ref_counts), alt_counts
        else:
            return None


    @classmethod
    def _get_variant_positions_from_vcf(cls, vcf_file):
        if not os.path.exists(vcf_file):
            return []
        f = pyfastaq.utils.open_file_read(vcf_file)
        positions = [l.rstrip().split('\t')[0:2] for l in f if not l.startswith('#')]
        positions = [(t[0], int(t[1]) - 1) for t in positions]
        pyfastaq.utils.close(f)
        return positions


    @staticmethod
    def _get_variants(vcf_file, read_depths_file, positions=None):
        if positions is None:
            positions = SamtoolsVariants._get_variant_positions_from_vcf(vcf_file)
        variants = {}
        if len(positions) == 0:
            return variants
        if not (os.path.exists(vcf_file) and os.path.exists(read_depths_file)):
            return variants
        for t in positions:
            name, pos = t[0], t[1]
            depths = SamtoolsVariants._get_read_depths(read_depths_file, name, pos)
            if depths is None:
                continue
            if name not in variants:
                variants[name] = {}
            variants[name][t[1]] = depths
        return variants


    @staticmethod
    def total_depth_per_contig(read_depths_file):
        f = pyfastaq.utils.open_file_read(read_depths_file)
        depths = {}
        for line in f:
            try:
                name, pos, base, var, depth, depth2 = line.rstrip().split('\t')
                depth = int(depth)
            except:
                pyfastaq.utils.close(f)
                raise Error('Error getting read depth from he following line of file ' + read_depths_file + ':\n' + line)

            depths[name] = depths.get(name, 0) + depth

        pyfastaq.utils.close(f)
        return depths


    @staticmethod
    def variants_in_coords(nucmer_matches, vcf_file):
        '''nucmer_matches = made by assembly_compare.assembly_match_coords().
           Returns number of variants that lie in nucmer_matches'''
        found_variants = {}
        f = pyfastaq.utils.open_file_read(vcf_file)
        for line in f:
            if line.startswith('#'):
                continue

            data = line.rstrip().split('\t')
            scaff = data[0]

            if scaff in nucmer_matches:
                position = int(data[1]) - 1
                i = pyfastaq.intervals.Interval(position, position)
                intersects = len([x for x in nucmer_matches[scaff] if x.intersects(i)]) > 0
                if intersects:
                    if scaff not in found_variants:
                        found_variants[scaff] = set()
                    found_variants[scaff].add(position)

        pyfastaq.utils.close(f)
        return found_variants


    def get_depths_at_position(self, seq_name, position):
        d = self._get_variants(self.vcf_file, self.read_depths_file, [(seq_name, position)])
        if seq_name in d and position in d[seq_name]:
            return d[seq_name][position]
        else:
            return 'ND', 'ND', 'ND', 'ND'


    def run(self):
        self._make_vcf_and_read_depths_files()
        # This is to make this object picklable, to keep multithreading happy
        self.log_fh = None
