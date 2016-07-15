class Error (Exception): pass

class Histogram:
    def __init__(self, bin_width=10):
        self.bin_width = bin_width
        self.bins = {}


    def __len__(self):
        if len(self.bins) == 0:
            return 0

        return sum(list(self.bins.values()))


    def __eq__(self, other):
       return type(other) is type(self) and self.__dict__ == other.__dict__


    def _to_bin(self, n):
        return self.bin_width * (n // self.bin_width)


    def add(self, n, count=1):
        b = self._to_bin(n)
        self.bins[b] = self.bins.get(b, 0) + count


    def stats(self):
        if len(self.bins) == 0:
            return None

        total = sum([self.bins[x] for x in self.bins])

        pc5 = None
        pc95 = None
        median = None
        sspace_sd = None
        i = 0

        for x in sorted(self.bins):
            i += self.bins[x]
            if i >= 0.05 * total and pc5 is None:
                pc5 = x + 0.5 * self.bin_width
            if i >= 0.5 * total and median is None:
                median = x + 0.5 * self.bin_width
            if i >= 0.95 * total and pc95 is None:
                pc95 = x + 0.5 + self.bin_width

        if None not in [pc5, pc95, median]:
            sspace_sd = round(max(pc95 - median, median - pc5) / median, 2)

        return (pc5, median, pc95, sspace_sd)

