"""Do true relatives stand out from a query's own background? For each query:
the top kinship values with their source, plus the 50th-highest value (the
query's breed-mix background). Truth: sib quartet DOGS-Gen-112..115, greyhound
littermates 79/80 and 82/83, first-degree DOGS-Gen-40/44. Artifact cases:
F1 crosses of inbred breeds (Cavalier/Bernese/Bichon x Poodle)."""
import gzip, json, sys
import numpy as np
D = '/u/project/pellegrini/gkislik/dogs'; REF = f'{D}/relatives_ref'
R = np.load(f'{REF}/geno.npy'); meta = json.load(gzip.open(f'{REF}/meta.json.gz', 'rt'))['samples']
ids = {s['id'].lower(): j for j, s in enumerate(meta)}
def king(q):
    use = q >= 0; Rs, qs = R[use], q[use]; valid = Rs >= 0
    het, hA, hB = qs == 1, qs == 0, qs == 2
    Nhh = (Rs[het] == 1).sum(0); Nopp = (Rs[hA] == 2).sum(0) + (Rs[hB] == 0).sum(0)
    den = valid[het].sum(0) + (Rs == 1).sum(0)
    return np.where(den > 0, (Nhh - 2.0 * Nopp) / np.maximum(den, 1), -1.0)
for name in sys.argv[1:]:
    j = ids.get(name.lower())
    if j is None: print(name, 'not in reference'); continue
    phi = king(R[:, j]); phi[j] = -1; phi[phi >= 0.45] = -1
    o = np.argsort(phi)[::-1]
    bg50, bg20 = phi[o[49]], phi[o[19]]
    top = [(round(float(phi[k]), 3), meta[k]['source'][:5], meta[k]['label'][:22], meta[k]['id'] if meta[k]['source'] != 'dap' else '') for k in o[:6]]
    print(f'{name:22s} {meta[j]["label"][:24]:24s} 20th {bg20:.3f} 50th {bg50:.3f} | ' + ' ; '.join(f'{p} {s} {l} {i}' for p, s, l, i in top))
