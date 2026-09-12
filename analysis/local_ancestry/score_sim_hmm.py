#!/usr/bin/env python3
# Score lai_hmm.py on the simulated admixed dogs (sim/<tag>.vcf.gz + .truth.tsv.gz
# from sim_admixed.py). The prior is the TRUE proportions (the lasso's job);
# the HMM sees only the genotypes.
#   python3 score_sim_hmm.py <sim_dir> [--temper T] [--gen auto|N] [--eps E] [--only substring]
import gzip, json, os, subprocess, sys, tempfile
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
args = sys.argv[1:]; simdir = args[0]
def opt(n, d): return args[args.index(n) + 1] if n in args else d
temper, gen, eps, only = opt('--temper', '10'), opt('--gen', 'auto'), opt('--eps', '0.05'), opt('--only', '')
panel = opt('--panel', '/Users/matteopellegrini/Downloads/dogs/breed_panel')
sites = []
with open(f'{panel}/sites.tsv') as f:
    next(f)
    for l in f: c, p = l.split()[:2]; sites.append((c, p))
site_idx = {s: i for i, s in enumerate(sites)}
results = []
for fn in sorted(os.listdir(simdir)):
    if not fn.endswith('.truth.tsv.gz') or only not in fn: continue
    tag = fn[:-len('.truth.tsv.gz')]
    truth = [l.rstrip('\n').split('\t') for l in gzip.open(f'{simdir}/{fn}', 'rt')]
    # truth rows are in the reference VCF site order == panel order (both from sites.tsv); assert lengths
    assert len(truth) == len(sites), (len(truth), len(sites))
    # true proportions -> a breed_result-shaped prior file
    cnt = defaultdict(int)
    for a, b in truth: cnt[a] += 1; cnt[b] += 1
    tot = sum(cnt.values())
    bj = {'breed_composition': [{'breed': b, 'proportion': n / tot} for b, n in cnt.items()]}
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as tf: json.dump(bj, tf); bjp = tf.name
    out = f'{simdir}/{tag}.hmm'
    r = subprocess.run([sys.executable, f'{HERE}/lai_hmm.py', f'{simdir}/{tag}.vcf.gz', bjp, out, '--temper', temper, '--gen', gen, '--eps', eps, '--floor', '0.001'], capture_output=True, text=True)
    if not os.path.exists(out + '.lai.json'): print(tag, 'FAILED', r.stderr[-300:]); continue
    res = json.load(open(out + '.lai.json'))
    # per-site unordered-pair accuracy from segments
    call = [None] * len(sites)
    for s in res['segments']:
        pass
    # rebuild per-site calls: segments are contiguous in site order per chromosome
    i = 0; segs = res['segments']; si = 0
    for m, (c, p) in enumerate(sites):
        while si < len(segs) - 1 and not (segs[si]['chrom'] == c and segs[si]['start'] <= int(p) <= segs[si]['end']): si += 1
        call[m] = tuple(sorted(segs[si]['anc']))
    ok = sum(1 for m in range(len(sites)) if call[m] == tuple(sorted(truth[m])))
    acc = ok / len(sites)
    tp = {b: n / tot for b, n in cnt.items()}
    mae = sum(abs(res['proportions'].get(b, 0) - tp.get(b, 0)) for b in set(tp) | set(res['proportions'])) / len(set(tp) | set(res['proportions']))
    true_segs = sum(1 for h in (0, 1) for m in range(len(sites)) if m == 0 or truth[m][h] != truth[m-1][h] or sites[m][0] != sites[m-1][0])
    results.append(dict(tag=tag, acc=acc, mae=mae, segs=res['n_segments'], true_segs=true_segs, gen=res['gen']))
    print(f'{tag:45s} site-acc {100*acc:5.1f}%  prop-MAE {100*mae:4.1f}pp  segments {res["n_segments"]:4d} (truth {true_segs})  gen {res["gen"]}')
print('\n=== summary (mean over replicates; temper', temper, 'gen', gen, 'eps', eps, ') ===')
groups = defaultdict(list)
for r in results:
    sc, _, noise = r['tag'].rsplit('.', 2); groups[(sc, noise)].append(r)
for (sc, noise), rs in groups.items():
    print(f'{sc:26s} {noise:14s} site-acc {100*sum(r["acc"] for r in rs)/len(rs):5.1f}%  prop-MAE {100*sum(r["mae"] for r in rs)/len(rs):4.1f}pp  segs {sum(r["segs"] for r in rs)/len(rs):.0f} vs truth {sum(r["true_segs"] for r in rs)/len(rs):.0f}')
