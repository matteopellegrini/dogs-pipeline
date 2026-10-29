#!/usr/bin/env python3
"""
Breed panel v2.1 source: SRA top-up genomes (GLIMPSE-imputed like customer dogs) -> npz on the Parker grid.

    python3 sra_dosage.py DS_DIR SHEET OUT.npz [--drop RUN ...]

DS_DIR holds <RUN>.tsv.gz (CHROM POS REF ALT DS at the panel sites, cluster/sra-ds-extract.sh); SHEET is
sample_sheet.sra.tsv (label = first word of the notes column). Hard-calls round(DS) and counts the Parker a1
allele (-1 where the site is missing or neither allele is a1), the same encoding as the other sources. Labels go
through SRA_ALIAS so they land on the panel's existing population names. Then removes duplicate animals
(IBS >= 0.97) against Parker, Dog10K and within the SRA set — several deposits re-release Dog10K/NHGRI genomes —
and writes the duplicates to sra_duplicates.tsv next to OUT.
"""
import glob, gzip, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from merged_panel import read_parker, read_dog10k

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
S = os.environ.get('SCRATCH', '/private/tmp/claude-501/-Users-matteopellegrini-Downloads-dogs/0fb26461-2afd-4142-936c-a58295d35959/scratchpad')
SRA_ALIAS = {'ANATOLIAN_SHEPHERD_DOG': 'ANATOLIAN_SHEPHERD', 'AUSTRALIAN_CATTLE_DOG': 'AUSTRALIAN_CATTLE',
             'CHINESE_SHAR_PEI': 'SHAR_PEI', 'GERMAN_SHEPHERD_DOG': 'GERMAN_SHEPHERD',
             'CATAHOULA_LEOPARD_DOG': 'CATAHOULA_LEOPARD', 'BOURBONNAIS_POINTING_DOG': 'BOURBONNAIS_POINTING'}
ds_dir, sheet, out = sys.argv[1], sys.argv[2], sys.argv[3]
drop = set(sys.argv[sys.argv.index('--drop') + 1:]) if '--drop' in sys.argv else set()

P, p_ids, p_truth, keys, a1 = read_parker(f'{ROOT}/COSMO/analysis/cosmo_parker_full', f'{ROOT}/COSMO/analysis/parker_clust.txt')
key_index = {k: i for i, k in enumerate(keys)}
label = {}
for l in open(sheet):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'sample_id' or len(f) < 9: continue
    run = os.path.basename(f[1]); lab = f[8].split()[0] if f[8].strip() else ''
    label[run] = SRA_ALIAS.get(lab, lab)
runs = sorted(r for r in (os.path.basename(p).split('.')[0] for p in glob.glob(f'{ds_dir}/*.tsv.gz')) if r in label and r not in drop)
dos = np.full((len(keys), len(runs)), -1, dtype=np.int8)
for j, r in enumerate(runs):
    for l in gzip.open(f'{ds_dir}/{r}.tsv.gz', 'rt'):
        c, pos, ref, alt, d = l.rstrip('\n').split('\t')[:5]
        i = key_index.get((c, pos))
        if i is None or ',' in alt: continue
        try: g = int(round(float(d)))
        except ValueError: continue
        if ref == a1[i]: dos[i, j] = 2 - g
        elif alt == a1[i]: dos[i, j] = g
print(f'{len(runs)} genomes loaded; call rate median {np.median((dos >= 0).mean(0)):.3f}', file=sys.stderr)

def ibs_matrix(A, B):
    okA = (A >= 0).mean(1) >= 0.95; okB = (B >= 0).mean(1) >= 0.95
    use = np.where(okA & okB)[0]
    A = A[use].astype(np.float32); B = B[use].astype(np.float32)
    mA = (A >= 0); mB = (B >= 0); A[~mA] = 0; B[~mB] = 0
    n_both = mA.T.astype(np.float32) @ mB.astype(np.float32)
    same = sum(((A == g) & mA).astype(np.float32).T @ ((B == g) & mB).astype(np.float32) for g in (0, 1, 2))
    off2 = ((A == 0) & mA).astype(np.float32).T @ ((B == 2) & mB).astype(np.float32) + ((A == 2) & mA).astype(np.float32).T @ ((B == 0) & mB).astype(np.float32)
    return (2 * same + (n_both - same - off2)) / (2 * np.maximum(n_both, 1))
D, d_ids, d_truth = read_dog10k(f'{S}/dog10k_gt.txt', f'{S}/dog10k_samples.txt', key_index, a1)
dup, rows = set(), []
for name, ibs, b_ids, b_lab in (('Parker', ibs_matrix(dos, P), p_ids, lambda i: p_truth.get(i, '')),
                                ('Dog10K', ibs_matrix(dos, D), d_ids, lambda i: d_truth.get(i, ''))):
    best = ibs.max(1); print(f'{name}: max-IBS quantiles 50/90/99% {np.quantile(best, [.5, .9, .99]).round(3)}', file=sys.stderr)
    for i, j in np.argwhere(ibs >= 0.97):
        dup.add(runs[i]); rows.append((name, runs[i], label[runs[i]], b_ids[j], b_lab(b_ids[j]), ibs[i, j]))
ss = ibs_matrix(dos, dos); np.fill_diagonal(ss, 0)
for i, j in np.argwhere(np.triu(ss) >= 0.97):
    if runs[i] not in dup and runs[j] not in dup:
        dup.add(runs[j]); rows.append(('SRA internal', runs[j], label[runs[j]], runs[i], label[runs[i]], ss[i, j]))
with open(os.path.join(os.path.dirname(os.path.abspath(out)), 'sra_duplicates.tsv'), 'w') as fh:
    fh.write('comparison\tsra_run\tsra_label\tother_id\tother_label\tibs\n')
    for r in rows: fh.write('\t'.join(map(str, r[:5])) + f'\t{r[5]:.4f}\n')
keep = [j for j, r in enumerate(runs) if r not in dup]
print(f'duplicates removed: {len(dup)}; kept {len(keep)}', file=sys.stderr)
np.savez_compressed(out, dos=dos[:, keep], ids=np.array([runs[j] for j in keep]), labels=np.array([label[runs[j]] for j in keep]))
