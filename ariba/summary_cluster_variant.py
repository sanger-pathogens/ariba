class Error (Exception): pass

class SummaryClusterVariant:
    def __init__(self, data_dict):
        self.known = None
        self.var_group = None
        self.coding = None
        self.var_string = None
        self.het_percent = None
        self._get_nonsynon_variant_data(data_dict)


    def __eq__(self, other):
       return type(other) is type(self) and self.__dict__ == other.__dict__


    def __hash__(self):
        return hash(tuple([self.__dict__[x] for x in sorted(self.__dict__.keys())]))


    def __str__(self):
        return ', '.join((str(self.known), self.var_group, str(self.coding), self.var_string, str(self.het_percent)))


    @classmethod
    def _has_nonsynonymous(cls, data_dict):
        return data_dict['ref_ctg_effect'] != 'SYN' and \
          (
              data_dict['has_known_var'] == '1' or \
              (data_dict['known_var'] != '1' and (data_dict['ref_ctg_change'] != '.' or data_dict['ref_ctg_effect'] != '.'))
          )


    @classmethod
    def _depths_look_het(cls, depths):
        sorted_depths = sorted(depths)
        min_var_read_depth = 5
        min_second_var_read_depth = 2
        max_allele_freq = 0.90
        total_depth = sum(depths)
        second_depth_ok = len(depths) == 1 or (len(depths) > 1 and sorted_depths[-2] >= min_second_var_read_depth)
        max_depth_ok = total_depth >= min_var_read_depth and sorted_depths[-1] / total_depth <= max_allele_freq
        return depths[0] < sorted_depths[-1] or (second_depth_ok and max_depth_ok)


    @classmethod
    def _get_is_het_and_percent(cls, data_dict):
        if data_dict['gene'] == '1' or not (data_dict['known_var'] == '1' or data_dict['ref_ctg_effect'] == 'SNP' or data_dict['var_type'] == 'HET') or data_dict['smtls_nts'] == '.' or ';' in data_dict['smtls_nts'] or data_dict['smtls_nts_depth'] == 'ND':
            return False, None
        else:
            nucleotides = data_dict['smtls_nts'].split(',')
            depths = data_dict['smtls_nts_depth'].split(',')

            if len(nucleotides) != len(depths):
                raise Error('Mismatch in number of inferred nucleotides from ctg_nt, smtls_nts, smtls_nts_depth columns. Cannot continue\n' + str(data_dict))

            try:
                depths = [int(x) for x in depths]
                nuc_to_depth = dict(zip(nucleotides, depths))

                if data_dict['ref_ctg_change'] != '.':
                    var_nucleotide = data_dict['ref_ctg_change'][-1]
                elif data_dict['var_type'] == 'HET':
                    var_nucleotide = '.'
                    best_depth = -1
                    for nuc in nuc_to_depth:
                        if nuc == data_dict['ctg_nt']:
                            continue
                        elif nuc_to_depth[nuc] > best_depth:
                            var_nucleotide = nuc
                            best_depth = nuc_to_depth[nuc]
                elif data_dict['known_var_change'] != '.':
                    var_nucleotide = data_dict['known_var_change'][-1]
                else:
                    return False, None

                if var_nucleotide == '.':
                    return False, None
                total_depth = sum(depths)
                if max([len(x) for x in nucleotides]) == 1:
                    var_depth = nuc_to_depth.get(var_nucleotide, 0)
                else:
                    var_depth = sum([nuc_to_depth[x] for x in nuc_to_depth if x[0] == var_nucleotide])

                is_het = SummaryClusterVariant._depths_look_het(depths)

                return is_het, round(100 * var_depth / total_depth, 1)
            except:
                return False, None


    def _get_nonsynon_variant_data(self, data_dict):
        self.known = data_dict['known_var'] == '1'
        self.coding = data_dict['gene'] == '1'
        self.var_group = data_dict['var_group']

        if data_dict['known_var'] == '1' and data_dict['known_var_change'] != '.':
            self.var_string = data_dict['known_var_change']
        elif data_dict['ref_ctg_change'] != '.':
            self.var_string = data_dict['ref_ctg_change']
        else:
            self.var_string = data_dict['ref_ctg_effect']

        self.is_het, self.het_percent = SummaryClusterVariant._get_is_het_and_percent(data_dict)

        if not SummaryClusterVariant._has_nonsynonymous(data_dict):
            self.has_nonsynon = False
            return

        self.has_nonsynon = True

        if data_dict['known_var_change'] == data_dict['ref_ctg_change'] == '.' == data_dict['ref_ctg_effect']:
            raise Error('Unexpected data in ariba summary... \n' + str(data_dict) + '\n... known_var_change, ref_ctg_change, ref_ctg_effect all equal to ".", but has a non synonymous change. Something is inconsistent. Cannot continue')
        elif '.' not in [data_dict['known_var_change'], data_dict['ref_ctg_change']] and \
          data_dict['known_var_change'] != data_dict['ref_ctg_change']:
            raise Error('Unexpected data in ariba summary... \n' + str(data_dict) + '\n... known_var_change != ref_ctg_change. Cannot continue')


