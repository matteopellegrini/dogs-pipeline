#!/usr/bin/env python3
"""Chunked PRS scorer for the Darwin's Ark bed (too big to hold in RAM).

    score_bed.py BFILE_PREFIX WEIGHTS.tsv OUT.tsv

WEIGHTS.tsv: snpid<TAB>beta (subset of bim SNPs). Streams the bed in blocks,
accumulates sum(beta * dosage) per dog with per-SNP mean imputation of
missing genotypes (dosage of A1).
"""
import sys
import numpy as np

def main():
    bfile, wfile, ofile = sys.argv[1], sys.argv[2], sys.argv[3]
    fam = [l.split() for l in open(bfile+'.fam')]
    n = len(fam)
    w = {}
    for l in open(wfile):
        s, b = l.split()
        w[s] = float(b)
    # row index of each weighted SNP
    idx, betas = [], []
    for i, l in enumerate(open(bfile+'.bim')):
        s = l.split()[1]
        if s in w:
            idx.append(i); betas.append(w[s])
    idx = np.array(idx); betas = np.array(betas, dtype=np.float64)
    bpf = (n+3)//4
    scores = np.zeros(n)
    LUT = np.zeros((256,4), dtype=np.float32)   # byte -> 4 dosages, -1=missing
    for byte in range(256):
        for k in range(4):
            c = (byte >> (2*k)) & 0b11
            LUT[byte,k] = 2 if c==0 else (1 if c==2 else (0 if c==3 else -1))
    with open(bfile+'.bed','rb') as fh:
        assert fh.read(3) == b'\x6c\x1b\x01'
        order = np.argsort(idx)
        pos_sorted = idx[order]; b_sorted = betas[order]
        CH = 20000
        for start in range(0, len(pos_sorted), CH):
            rows = pos_sorted[start:start+CH]
            bts = b_sorted[start:start+CH]
            fh.seek(3 + int(rows[0])*bpf)
            span = int(rows[-1]) - int(rows[0]) + 1
            buf = np.frombuffer(fh.read(span*bpf), dtype=np.uint8).reshape(span, bpf)
            sel = buf[rows - rows[0]]
            dos = LUT[sel].reshape(len(rows), -1)[:, :n]
            miss = dos < 0
            if miss.any():
                dsum = np.where(miss, 0, dos).sum(axis=1)
                cnt = (~miss).sum(axis=1).clip(min=1)
                mu = (dsum/cnt)[:,None]
                dos = np.where(miss, mu, dos)
            scores += bts @ dos
    with open(ofile,'w') as out:
        for f, s in zip(fam, scores):
            out.write(f'{f[0]}\t{s:.6f}\n')
    print(f'{len(idx)}/{len(w)} weighted SNPs found in bim; scored {n} dogs')

if __name__ == '__main__':
    main()
