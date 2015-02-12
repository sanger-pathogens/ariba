import pyfastaq
from ariba import link

class Error (Exception): pass

class Graph:
    def __init__(self, ref_lengths):
        self.links = {}
        self.partial_links = {}
        self.ref_lengths = ref_lengths


    def update_from_sam(self, sam, sam_reader):
        '''Updates graph info from a pysam.AlignedSegment object'''
        if sam.is_unmapped \
          or sam.mate_is_unmapped \
          or (sam.reference_id == sam.next_reference_id):
            return

        new_link = link.Link(sam, sam_reader, self.ref_lengths)
        read_name = sam.query_name

        if read_name in self.partial_links:
            new_link.merge(self.partial_links[read_name])
            del self.partial_links[read_name]
            key = tuple(sorted((new_link.refnames[0], new_link.refnames[1])))
            if key not in self.links:
                self.links[key] = []
            new_link.sort()
            self.links[key].append(new_link)
        else:
            self.partial_links[read_name] = new_link


    def _make_graph(self, max_insert):
        '''helper function to construct graph from current state of object'''
        if len(self.partial_links) != 0:
            raise Error('Error in _make_graph(). Cannot continue because there are partial links')
        
        self.contig_links = {}
        for key in self.links:
            for l in self.links[key]:
                insert_size = l.insert_size()
                if insert_size <= max_insert:
                    if key not in self.contig_links:
                        self.contig_links[key] = {}
                    dirs = ''.join(l.dirs)
                    self.contig_links[key][dirs] = self.contig_links[key].get(dirs, 0) + 1


    def _remove_low_cov_links(self, min_coverage):
        keys_to_delete = set()

        for key in self.contig_links:
            d = self.contig_links[key]
            self.contig_links[key] = {x:d[x] for x in d if d[x] >= min_coverage}
            if len(self.contig_links[key]) == 0:
                keys_to_delete.add(key)

        for key in keys_to_delete:
            del self.contig_links[key]


    def _contig_graph_is_consistent(self):
        for key in self.contig_links:
            if len(self.contig_links[key]) > 1:
                return False

        all_names_in_keys = set()
        for key in self.contig_links:
            all_names_in_keys.add(key[0])
            all_names_in_keys.add(key[1])

        name_to_links = {name: set() for name in all_names_in_keys}

        for key in self.contig_links:
            assert len(self.contig_links[key]) == 1
            directions = list(self.contig_links[key])[0]
            name1, name2 = key
            dir1, dir2 = directions[0], directions[1]
            for n, d in [(name1, dir1), (name2, dir2)]:
                if d in name_to_links[n]:
                    return False
                else:
                    name_to_links[n].add(d)
        return True


    def is_consistent(self, min_coverage, max_insert):
        self._make_graph(max_insert)
        self._remove_low_cov_links(min_coverage)
        return self._contig_graph_is_consistent()


    def write_all_links_to_file(self, filename):
        f = pyfastaq.utils.open_file_write(filename)
        for key in sorted(self.links):
            to_write = list(self.links[key])
            to_write.sort()
            for l in to_write:
                print(l, file=f)
        pyfastaq.utils.close(f)

