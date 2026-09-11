# Two-step local ancestry for one dog (pilot recipe, 2026-09-11):
#   lasso proportions -> breeds >= FLOOR -> reference restricted to them ->
#   FLARE with those proportions as prior, EM off, gen=GEN -> genome-wide
#   proportions (.global.anc.gz) + painting (.anc.vcf.gz).
import gzip, json, subprocess, sys, time
q, breed_json, tag = sys.argv[1:4]
FLOOR = float(sys.argv[4]) if len(sys.argv) > 4 else 0.02
GEN = sys.argv[5] if len(sys.argv) > 5 else '3'
J = '/opt/homebrew/opt/openjdk/bin/java'
# lasso -> prior over ALL reference breeds (first-appearance order), then floor
rows = [l.rstrip('\n').split('\t') for l in open('ref_panel_map.tsv') if l.strip()]
b = json.load(open(breed_json)); prior = {}
for e in b['breed_composition']:
    for c in (e.get('components') or [{'code': e.get('breed'), 'proportion': e.get('proportion', 0)}]):
        prior[c.get('code') or e.get('breed')] = prior.get(c.get('code') or e.get('breed'), 0.0) + float(c.get('proportion', 0))
refb = {bb for _, bb in rows}
keep = {bb: v for bb, v in prior.items() if bb in refb and v >= FLOOR}
tot = sum(keep.values()); keep = {bb: v / tot for bb, v in keep.items()}
order = []; ks = set()
for s, bb in rows:
    if bb in keep: ks.add(s); (order.append(bb) if bb not in order else None)
sid = subprocess.run(f"gzip -dc {q} | grep -m1 '^#CHROM' | cut -f10", shell=True, capture_output=True, text=True).stdout.strip()
open(f'work_{tag}.map.tsv', 'w').write(''.join(f'{s}\t{bb}\n' for s, bb in rows if bb in keep))
open(f'work_{tag}.prior.tsv', 'w').write('SAMPLE ' + ' '.join(order) + '\n' + sid + ' ' + ' '.join(f'{keep[bb]:.6f}' for bb in order) + '\n')
with gzip.open('ref_panel.vcf.gz', 'rt') as f, gzip.open(f'work_{tag}.ref.vcf.gz', 'wt') as o:
    for l in f:
        if l.startswith('##'): o.write(l); continue
        p = l.rstrip('\n').split('\t')
        if l.startswith('#'): idx = [i for i in range(9, len(p)) if p[i] in ks]; o.write('\t'.join(p[:9] + [p[i] for i in idx]) + '\n'); continue
        o.write('\t'.join(p[:9] + [p[i] for i in idx]) + '\n')
def flare_run(tag, keep, order, ks):
    t = time.time()
    r = subprocess.run([J, '-Xmx8g', '-jar', 'flare.jar', f'ref=work_{tag}.ref.vcf.gz', f'ref-panel=work_{tag}.map.tsv', f'gt={q}',
                    'map=constant_1cM_per_Mb.map', f'out=out/{tag}', 'probs=true', f'gt-ancestries=work_{tag}.prior.tsv',
                    'em=false', f'gen={GEN}', 'nthreads=8', 'seed=1'], capture_output=True, text=True)
    gl = [l.split() for l in gzip.open(f'out/{tag}.global.anc.gz', 'rt')]
    return dict(zip(gl[0][1:], map(float, gl[1][1:]))), time.time() - t
props, secs = flare_run(tag, keep, order, ks)
# Sink guard: a small lasso component (<5%) whose painted share balloons past 3x
# its weight is absorbing segments from a related breed (ZAGAR at low depth);
# drop it, rebuild the reference/prior, and run once more.
sinks = [bb for bb, v in props.items() if keep.get(bb, 0) < 0.05 and v > 3 * keep.get(bb, 1e-9)]
if sinks:
    for bb in sinks: keep.pop(bb, None)
    tot = sum(keep.values()); keep = {bb: v / tot for bb, v in keep.items()}
    order = []; ks = set()
    for s_, bb in rows:
        if bb in keep: ks.add(s_); (order.append(bb) if bb not in order else None)
    open(f'work_{tag}.map.tsv', 'w').write(''.join(f'{s_}\t{bb}\n' for s_, bb in rows if bb in keep))
    open(f'work_{tag}.prior.tsv', 'w').write('SAMPLE ' + ' '.join(order) + '\n' + sid + ' ' + ' '.join(f'{keep[bb]:.6f}' for bb in order) + '\n')
    with gzip.open('ref_panel.vcf.gz', 'rt') as f, gzip.open(f'work_{tag}.ref.vcf.gz', 'wt') as o:
        for l in f:
            if l.startswith('##'): o.write(l); continue
            p = l.rstrip('\n').split('\t')
            if l.startswith('#'): idx = [i for i in range(9, len(p)) if p[i] in ks]; o.write('\t'.join(p[:9] + [p[i] for i in idx]) + '\n'); continue
            o.write('\t'.join(p[:9] + [p[i] for i in idx]) + '\n')
    props, secs = flare_run(tag, keep, order, ks)
    print(f'  [sink guard] pruned {sinks} and re-ran')
nseg = 0; prev = {}; chrom = None
for l in gzip.open(f'out/{tag}.anc.vcf.gz', 'rt'):
    if l.startswith('#'): continue
    p = l.split('\t', 10); fmt = p[8].split(':'); v = p[9].split(':')
    if p[0] != chrom: chrom = p[0]; prev = {}
    for h, a in ((1, v[fmt.index('AN1')]), (2, v[fmt.index('AN2')])):
        if a != prev.get(h): nseg += 1; prev[h] = a
print(f'{tag}: {len(order)} breeds / {len(idx)} ref dogs, FLARE {secs:.0f}s, {nseg} segments; '
      + ', '.join(f'{bb} {v:.3f} (lasso {keep[bb]:.3f})' for bb, v in sorted(props.items(), key=lambda kv: -kv[1])[:6]))
json.dump({'tag': tag, 'flare': props, 'lasso_floor': keep, 'segments': nseg, 'pruned_sinks': sinks}, open(f'out/{tag}.summary.json', 'w'))
