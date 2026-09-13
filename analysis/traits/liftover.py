#!/usr/bin/env python3
"""
Minimal UCSC chain-file liftover (no external deps) for single positions.

  from liftover import Lifter
  L = Lifter('darwins_ark/canFam3ToCanFam4.over.chain.gz')
  L.convert('chr32', 4509367)   # 1-based in -> (chrom, 1-based pos, strand) or None

  python3 liftover.py <chain.gz> chr32:4509367 chr17:19431807 ...
"""
import gzip, sys, bisect


class Lifter:
    def __init__(self, chain_path):
        self.blocks = {}   # tchrom -> list of (tstart, tend, qchrom, qstart, qstrand, qsize)  (0-based half-open on target)
        op = gzip.open if chain_path.endswith('.gz') else open
        with op(chain_path, 'rt') as f:
            hdr = None
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith('chain'):
                    p = line.split()
                    # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                    hdr = dict(tname=p[2], tstart=int(p[5]), qname=p[7], qsize=int(p[8]), qstrand=p[9], qstart=int(p[10]))
                    tpos, qpos = hdr['tstart'], hdr['qstart']
                    continue
                p = line.split()
                size = int(p[0])
                self.blocks.setdefault(hdr['tname'], []).append((tpos, tpos + size, hdr['qname'], qpos, hdr['qstrand'], hdr['qsize']))
                if len(p) == 3:
                    tpos += size + int(p[1]); qpos += size + int(p[2])
        for c in self.blocks:
            self.blocks[c].sort()
        self.starts = {c: [b[0] for b in v] for c, v in self.blocks.items()}

    def convert(self, chrom, pos1):
        """1-based position -> (qchrom, 1-based qpos, strand) or None"""
        t = pos1 - 1
        bl = self.blocks.get(chrom)
        if not bl: return None
        i = bisect.bisect_right(self.starts[chrom], t) - 1
        if i < 0: return None
        tstart, tend, qchrom, qstart, qstrand, qsize = bl[i]
        if not (tstart <= t < tend): return None
        off = t - tstart
        if qstrand == '+':
            return qchrom, qstart + off + 1, '+'
        # minus strand: query coordinates count from the reverse-complemented end
        q = qsize - (qstart + off) - 1
        return qchrom, q + 1, '-'


if __name__ == '__main__':
    L = Lifter(sys.argv[1])
    for a in sys.argv[2:]:
        c, p = a.split(':'); r = L.convert(c, int(p))
        print(a, '->', f'{r[0]}:{r[1]} ({r[2]})' if r else 'unmapped')
