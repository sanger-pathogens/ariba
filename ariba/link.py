class Error (Exception): pass

class Link:
    def __init__(self, sam, sam_reader, ref_lengths, s=None):
        # s meant for writing tests. Construct object from a string (same format as __str__())
        if s is not None:
            l = s.split('\t')
            assert len(l) == 8
            self.refnames = [l[0], l[4]]
            self.lengths = [None if l[i] == '.' else int(l[i]) for i in [1,5]]
            self.dirs = [l[2], l[6]]
            self.pos = [None if l[i] == '.' else int(l[i]) for i in [3,7]]
            return

        if sam.is_unmapped or sam.mate_is_unmapped or sam.reference_id == sam.next_reference_id:
            raise Error('Cannot create link from sam:\n' + str(sam))

        self.refnames = [sam_reader.getrname(sam.reference_id), sam_reader.getrname(sam.next_reference_id)]
        self.dirs = ['L' if sam.is_reverse else 'R', 'L' if sam.mate_is_reverse else 'R']
        position = sam.reference_end - 1 if sam.is_reverse else sam.reference_start
        self.pos = [position, None]
        self.lengths = [ref_lengths[x] for x in self.refnames]
        
        if not sam.is_read1:
            self._swap()


    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__


    def __lt__(self, other):
        if self.refnames[0] < other.refnames[0]:
            return True
        elif self.refnames[0] > other.refnames[0]:
            return False
        elif self.pos[0] < other.pos[0]:
            return True
        elif self.pos[0] == other.pos[0]:
            return self.pos[1] < other.pos[1]
        else:
            return False


    def __str__(self):
        return '\t'.join([
           self.refnames[0], 
           str(self.lengths[0]),
           self.dirs[0],
           str(self.pos[0]) if self.pos[0] is not None else '.',
           self.refnames[1], 
           str(self.lengths[1]),
           self.dirs[1],
           str(self.pos[1]) if self.pos[1] is not None else '.'
        ])


    def _swap(self):
        for l in [self.refnames, self.dirs, self.pos, self.lengths]:
            l.reverse()


             

    def sort(self):
        if self.refnames[1] < self.refnames[0]:
            self._swap()


    def _distance_to_contig_end(self, one_or_two):
        one_or_two -= 1
        if self.dirs[one_or_two] == 'L' and self.pos[one_or_two] is not None:
            return self.pos[one_or_two]
        elif self.dirs[one_or_two] == 'R' and self.pos[one_or_two] is not None:
            return self.lengths[one_or_two] - self.pos[one_or_two] - 1
        else:
            raise Error('Error in _distance_to_contig_end(' + str(one_or_two + 1) + ') for link:\n' + str(self))
        

    def merge(self, other):
        '''Merge another link into this one. Expected that each link was created from each mate from a pair. We only know both distances to contig ends when we have read info from both mappings in a BAM file. All other info should be the same.'''
        assert self.refnames == other.refnames
        assert self.dirs == other.dirs
        assert self.lengths == other.lengths

        for i in range(2):
            if self.pos[i] is None:
                if other.pos[i] is None:
                   raise Error('Error merging these two links:\n' + str(self) + '\n' + str(other))
                self.pos[i] = other.pos[i]
            else:
                if other.pos[i] is not None:
                   raise Error('Error merging these two links:\n' + str(self) + '\n' + str(other))


    def _distance_to_contig_ends(self):
        return self._distance_to_contig_end(1), self._distance_to_contig_end(2)


    def insert_size(self):
        '''Returns insert size, defined as distance from outer edges of reads (and assumes gap length of zero)'''
        try:
            distances = self._distance_to_contig_ends()
        except:
            raise Error('Error getting insert size from Link:\n' + str(self))

        return sum(distances)
