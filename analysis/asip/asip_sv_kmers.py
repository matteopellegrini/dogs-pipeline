"""Count reads carrying each allele's breakpoint k-mers at the ASIP structural
differences between canFam3 (Boxer, dominant yellow VP1-HCP1) and canFam4
(German Shepherd). Events come from the minimap2 cs string of cf3 vs cf4.
usage: asip_sv_kmers.py cf3.fa cf4.fa cf3_vs_cf4.paf bam [bam ...]"""
import re, sys, pysam
K = 31
def fa(p): return ''.join(l.strip() for l in open(p) if not l.startswith('>')).upper()
cf3, cf4 = fa(sys.argv[1]), fa(sys.argv[2])
CF4_OFF, CF3_OFF = 23830000, 23290000          # 0-based genomic start of each slice
rc = lambda s: s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
events = []
for l in open(sys.argv[3]):
    f = l.split('\t'); cs = [x for x in f if x.startswith('cs:Z:')][0][5:]
    t, q = int(f[7]), int(f[2])
    for op in re.findall(r'(:\d+|\*[a-z]{2}|[+-][a-z]+)', cs):
        if op[0] == ':': n = int(op[1:]); t += n; q += n
        elif op[0] == '*': t += 1; q += 1
        else:
            n = len(op) - 1
            if n >= 15:
                if op[0] == '+':   # sequence in cf3 only
                    alt = [cf3[q - K//2:q + K//2 + 1], cf3[q + n - K//2:q + n + K//2 + 1]]
                    ref = [cf4[t - K//2:t + K//2 + 1]]
                    events.append(dict(kind='cf3_only', len=n, cf4=t + CF4_OFF, ref=ref, alt=alt))
                else:              # sequence in cf4 only
                    ref = [cf4[t - K//2:t + K//2 + 1], cf4[t + n - K//2:t + n + K//2 + 1]]
                    alt = [cf3[q - K//2:q + K//2 + 1]]
                    events.append(dict(kind='cf4_only', len=n, cf4=t + CF4_OFF, ref=ref, alt=alt))
            if op[0] == '+': q += n
            else: t += n
for bam in sys.argv[4:]:
    b = pysam.AlignmentFile(bam)
    print('#', bam)
    for e in events:
        kr = set(e['ref']) | {rc(x) for x in e['ref']}
        ka = set(e['alt']) | {rc(x) for x in e['alt']}
        nr = na = 0
        for r in b.fetch('chr24', e['cf4'] - 600, e['cf4'] + e['len'] + 600):
            if r.is_secondary or r.is_supplementary or r.is_duplicate: continue
            s = r.query_sequence or ''
            if any(k in s for k in ka): na += 1
            elif any(k in s for k in kr): nr += 1
        print(f"{e['kind']:9s} len {e['len']:5d} cf4 {e['cf4']:,}  ref-junction reads {nr:3d}  alt-junction reads {na:3d}")
