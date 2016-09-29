import os
import sys
import pysam
import pyfastaq
import vcfcall_ariba

class Error (Exception): pass


class SamtoolsVariants:
    def __init__(self,
      ref_fa,
      bam,
      outprefix,
      log_fh=sys.stdout,
      min_var_read_depth=5,
      min_second_var_read_depth=2,
      max_allele_freq=0.90
    ):
        self.ref_fa = os.path.abspath(ref_fa)
        self.bam = os.path.abspath(bam)
        self.outprefix = os.path.abspath(outprefix)
        self.log_fh = log_fh
        self.min_var_read_depth = min_var_read_depth
        self.min_second_var_read_depth = min_second_var_read_depth
        self.max_allele_freq = max_allele_freq

        self.vcf_file = self.outprefix + '.vcf'
        self.read_depths_file = self.outprefix + '.read_depths.gz'
        self.contig_depths_file = self.outprefix + '.contig_depths'


    def _make_vcf_and_read_depths_files(self):
        if not os.path.exists(self.ref_fa + '.fai'):
            pysam.faidx(self.ref_fa)

        tmp_vcf = self.vcf_file + '.tmp'
        with open(tmp_vcf, 'w') as f:
            print(pysam.mpileup(
                '-t', 'INFO/AD',
                '-L', '99999999',
                '-A',
                '-f', self.ref_fa,
                '-u',
                '-v',
                self.bam,
            ), end='', file=f)

        got = vcfcall_ariba.vcfcall_ariba(tmp_vcf, self.outprefix, self.min_var_read_depth, self.min_second_var_read_depth, self.max_allele_freq)
        if got != 0:
            raise Error('Error parsing vcf file. Cannot contine')

        pysam.tabix_compress(self.outprefix + '.read_depths', self.read_depths_file)
        pysam.tabix_index(self.read_depths_file, seq_col=0, start_col=1, end_col=1)
        os.unlink(self.outprefix + '.read_depths')
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
            bases = ref_base if alt_base == '.' else ref_base + ',' + alt_base
            return bases, int(ref_counts), alt_counts
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
                name, depth = line.rstrip().split('\t')
                depth = int(depth)
            except:
                pyfastaq.utils.close(f)
                raise Error('Error getting read depth from he following line of file ' + read_depths_file + ':\n' + line)

            depths[name] = depth

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
            return 'ND', 'ND', 'ND'


    def run(self):
        self._make_vcf_and_read_depths_files()
        # This is to make this object picklable, to keep multithreading happy
        self.log_fh = None
