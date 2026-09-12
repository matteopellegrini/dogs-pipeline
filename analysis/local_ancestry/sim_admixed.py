#!/usr/bin/env python3
# Simulated-admixture validation of the two-step local-ancestry recipe.
#
# For each scenario: hold out 2 dogs per source breed from the phased reference,
# splice their haplotypes into admixed dogs (segment lengths ~ Exp(100/gen cM)
# on the constant 1 cM/Mb map; F1 = one pure haplotype per parent breed), write
# a phased query VCF with optional noise (phase switches every ~S Mb; genotype
# error rate e), run FLARE with the scenario breeds as reference (minus the
# held-out dogs) and the TRUE proportions as prior (the lasso's job), and score
# per-site ancestry (unordered haplotype pair) against the truth.
import gzip, json, os, random, subprocess, sys
from collections import defaultdict
J = '/opt/homebrew/opt/openjdk/bin/java'
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 7)

SCEN = {
  'F1_poodle_x_cocker':      dict(props={'STANDARD_POODLE': .5, 'COCKER_SPANIEL': .5}, gen=1),
  'doodle_3gen':             dict(props={'STANDARD_POODLE': .5, 'MINIATURE_POODLE': .25, 'LABRADOR_RETRIEVER': .125, 'COCKER_SPANIEL': .125}, gen=3),
  'close_relatives_GSD_WSS': dict(props={'GERMAN_SHEPHERD': .5, 'WHITE_SWISS_SHEPHERD': .5}, gen=2),
  'supermutt_6':             dict(props={'BOXER': 1/6, 'BEAGLE': 1/6, 'CHIHUAHUA': 1/6, 'GOLDEN_RETRIEVER': 1/6, 'SIBERIAN_HUSKY': 1/6, 'AMERICAN_PIT_BULL_TERRIER': 1/6}, gen=5),
}
NOISE = {'perfect': (None, 0.0), 'phase5Mb': (5e6, 0.0), 'phase5Mb_err5': (5e6, 0.05)}
N_REP = 3
HOLD = 2

rows = [l.rstrip('\n').split('\t') for l in open('ref_panel_map.tsv') if l.strip()]
by_breed = defaultdict(list)
for s, b in rows: by_breed[b].append(s)

# --- load reference: header + per-site records (kept as text lines for speed) ---
hdr = None; sites = []; recs = []
with gzip.open('ref_panel.vcf.gz', 'rt') as f:
    for l in f:
        if l.startswith('##'): continue
        if l.startswith('#'): hdr = l.rstrip('\n').split('\t'); continue
        recs.append(l.rstrip('\n'))
col = {s: i for i, s in enumerate(hdr)}
print(f'reference loaded: {len(recs)} sites, {len(hdr)-9} dogs')
# chromosome -> list of (site index, pos)
chrom_sites = defaultdict(list)
for i, r in enumerate(recs):
    c, p = r.split('\t', 2)[:2]; chrom_sites[c].append((i, int(p)))

def hap_alleles(sample, hap):
    """allele (0/1 char) of one haplotype of one reference dog at every site"""
    ci = col[sample]; out = []
    for r in recs:
        g = r.split('\t')[ci].split(':')[0]
        out.append(g[0] if hap == 0 else g[2])
    return out

def simulate(props, gen, held):
    """returns (hap1_alleles, hap2_alleles, truth1, truth2) over all sites"""
    breeds = list(props); w = [props[b] for b in breeds]
    src = {b: [(s, h) for s in held[b] for h in (0, 1)] for b in breeds}
    cache = {}
    def alle(s, h):
        if (s, h) not in cache: cache[(s, h)] = hap_alleles(s, h)
        return cache[(s, h)]
    out = []
    for hp in range(2):
        alle_out = [None] * len(recs); truth = [None] * len(recs)
        if gen == 1:   # F1: whole haplotype from one parent breed
            b = breeds[hp % len(breeds)]; s, h = random.choice(src[b]); a = alle(s, h)
            for i in range(len(recs)): alle_out[i] = a[i]; truth[i] = b
        else:
            for c, lst in chrom_sites.items():
                pos0 = lst[0][1]; end = lst[-1][1]; cur = pos0
                while cur <= end:
                    b = random.choices(breeds, w)[0]; s, h = random.choice(src[b]); a = alle(s, h)
                    seg_len = random.expovariate(gen / 100.0) * 1e6   # cM -> bp at 1 cM/Mb
                    nxt = cur + seg_len
                    for i, p in lst:
                        if cur <= p < nxt: alle_out[i] = a[i]; truth[i] = b
                    cur = nxt
        out.append((alle_out, truth))
    return out[0][0], out[1][0], out[0][1], out[1][1]

def write_query(path, sid, h1, h2, switch_every, err):
    """phased VCF; phase switches swap h1/h2 from a random point; genotype errors flip an allele"""
    with gzip.open(path, 'wt') as o:
        o.write('##fileformat=VCFv4.2\n')
        for c in chrom_sites: o.write(f'##contig=<ID={c}>\n')
        o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased">\n')
        o.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sid + '\n')
        for c, lst in chrom_sites.items():
            flipped = False; nxt = lst[0][1] + (random.expovariate(1 / switch_every) if switch_every else 1e18)
            for i, p in lst:
                if switch_every and p >= nxt: flipped = not flipped; nxt = p + random.expovariate(1 / switch_every)
                a, b = (h1[i], h2[i]) if not flipped else (h2[i], h1[i])
                if err and random.random() < err: a = '1' if a == '0' else '0'
                if err and random.random() < err: b = '1' if b == '0' else '0'
                f = recs[i].split('\t', 5)
                o.write(f'{f[0]}\t{f[1]}\t{f[2]}\t{f[3]}\t{f[4]}\t.\tPASS\t.\tGT\t{a}|{b}\n')

def write_ref(path, keep):
    ks = set(keep); idx = [i for i in range(9, len(hdr)) if hdr[i] in ks]
    with gzip.open(path, 'wt') as o:
        o.write('##fileformat=VCFv4.2\n')
        for c in chrom_sites: o.write(f'##contig=<ID={c}>\n')
        o.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased">\n')
        o.write('\t'.join(hdr[:9] + [hdr[i] for i in idx]) + '\n')
        for r in recs:
            f = r.split('\t'); o.write('\t'.join(f[:9] + [f[i] for i in idx]) + '\n')

def score(prefix, truth1, truth2):
    names = {}; ok = n = nseg = 0; prev = {}; chrom = None; props = defaultdict(int)
    for l in gzip.open(prefix + '.anc.vcf.gz', 'rt'):
        if l.startswith('##ANCESTRY'):
            for a in l.strip().split('=', 1)[1].strip('<>').split(','):
                k, v = a.strip().split('='); names[int(v)] = k
            continue
        if l.startswith('#'): continue
        p = l.rstrip('\n').split('\t'); fmt = p[8].split(':'); g = p[9].split(':')
        a1, a2 = names[int(g[fmt.index('AN1')])], names[int(g[fmt.index('AN2')])]
        i = site_index[(p[0], p[1])]
        n += 1; ok += sorted((a1, a2)) == sorted((truth1[i], truth2[i]))
        props[a1] += 1; props[a2] += 1
        if p[0] != chrom: chrom = p[0]; prev = {}
        for h, a in ((1, a1), (2, a2)):
            if a != prev.get(h): nseg += 1; prev[h] = a
    tot = sum(props.values()); est = {b: v / tot for b, v in props.items()}
    return ok / n, est, nseg

site_index = {}
for i, r in enumerate(recs):
    c, p = r.split('\t', 2)[:2]; site_index[(c, p)] = i

os.makedirs('sim', exist_ok=True)
results = []
for sname, sc in SCEN.items():
    props, gen = sc['props'], sc['gen']
    held = {b: random.sample(by_breed[b], HOLD) for b in props}
    heldset = {s for v in held.values() for s in v}
    ref_keep = [s for s, b in rows if b in props and s not in heldset]
    write_ref(f'sim/{sname}.ref.vcf.gz', ref_keep)
    open(f'sim/{sname}.map.tsv', 'w').write(''.join(f'{s}\t{b}\n' for s, b in rows if s in set(ref_keep)))
    order = []
    for s, b in rows:
        if s in set(ref_keep) and b not in order: order.append(b)
    # truth of the simulated dog's SEGMENTS is what the lasso would estimate; use the scenario proportions as prior
    for rep in range(N_REP):
        h1, h2, t1, t2 = simulate(props, gen, held)
        true_prop = defaultdict(int)
        for t in t1 + t2: true_prop[t] += 1
        tp = {b: v / (2 * len(recs)) for b, v in true_prop.items()}
        for nname, (sw, err) in NOISE.items():
            tag = f'{sname}.r{rep}.{nname}'; sid = f'SIM_{tag}'
            write_query(f'sim/{tag}.vcf.gz', sid, h1, h2, sw, err)
            with gzip.open(f'sim/{tag}.truth.tsv.gz', 'wt') as tf:
                for i in range(len(recs)): tf.write(f'{t1[i]}\t{t2[i]}\n')
            if os.environ.get('SKIP_FLARE'): continue
            open(f'sim/{tag}.prior.tsv', 'w').write('SAMPLE ' + ' '.join(order) + '\n' + sid + ' ' + ' '.join(f'{tp.get(b, 0.001):.6f}' for b in order) + '\n')
            r = subprocess.run([J, '-Xmx8g', '-jar', 'flare.jar', f'ref=sim/{sname}.ref.vcf.gz', f'ref-panel=sim/{sname}.map.tsv',
                                f'gt=sim/{tag}.vcf.gz', 'map=constant_1cM_per_Mb.map', f'out=sim/{tag}', 'probs=true',
                                f'gt-ancestries=sim/{tag}.prior.tsv', 'em=false', f'gen={gen}', 'min-mac=2', 'nthreads=8', 'seed=1'],
                               capture_output=True, text=True)
            if not os.path.exists(f'sim/{tag}.anc.vcf.gz') or os.path.getsize(f'sim/{tag}.anc.vcf.gz') == 0:
                print(f'{tag}: FLARE FAILED: {r.stderr[-300:]}'); continue
            acc, est, nseg = score(f'sim/{tag}', t1, t2)
            mae = sum(abs(est.get(b, 0) - tp.get(b, 0)) for b in set(est) | set(tp)) / len(set(est) | set(tp))
            true_segs = sum(1 for t in (t1, t2) for i in range(1, len(t)) if t[i] != t[i-1] or recs[i].split('\t',1)[0] != recs[i-1].split('\t',1)[0]) + 2 * len(chrom_sites)
            results.append(dict(scenario=sname, rep=rep, noise=nname, site_acc=acc, prop_mae=mae, segs=nseg, true_segs=true_segs))
            print(f'{tag:45s} site-acc {100*acc:5.1f}%  prop-MAE {100*mae:4.1f}pp  segments {nseg} (truth {true_segs})')
json.dump(results, open('sim/results.json', 'w'), indent=1)
print('\n=== summary (mean over replicates) ===')
for sname in SCEN:
    for nname in NOISE:
        rs = [r for r in results if r['scenario'] == sname and r['noise'] == nname]
        if rs: print(f'{sname:26s} {nname:14s} site-acc {100*sum(r["site_acc"] for r in rs)/len(rs):5.1f}%  prop-MAE {100*sum(r["prop_mae"] for r in rs)/len(rs):4.1f}pp  segs {sum(r["segs"] for r in rs)/len(rs):.0f} vs truth {sum(r["true_segs"] for r in rs)/len(rs):.0f}')
