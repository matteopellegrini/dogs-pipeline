#!/usr/bin/env python3
"""
DLA (dog MHC) class II diversity from Dog10K-imputed, phased genotypes.

Layer 1 of the MHC feature (2026-09-13): every processed dog has GLIMPSE2
phased genotypes at the 5,546 Dog10K panel sites in the class II window
chr12:2,300,000-2,660,000 (canFam4: BTNL2 .. DLA-DRA .. DLA-DRB1/DQA1/DQB1 ..
TAP2). From those we score, per dog:

  het_rate        fraction of window sites that are heterozygous (imputed GT)
  hap_distance    Hamming distance (fraction of sites) between the dog's two
                  phased haplotypes in the CORE block (DRB1-DQA1-DQB1)
  hap_groups      the two haplotypes assigned to haplotype GROUPS defined on
                  the 1,929 breed-labelled Dog10K panel dogs (greedy clustering
                  at <= CLUSTER_THR mismatches in the core block); a dog is
                  "homozygous" for the class II haplotype when both copies fall
                  in the same group
  percentiles     against the cohort and against panel dogs of the same breed

Inputs (built by cluster/mhc_extract.sh + the panel query):
  mhc/classII_sites.tsv          chrom pos ref alt (window sites, panel order)
  mhc/panel_classII_gt.tsv       chrom pos ref alt GT... (one column per panel sample)
  mhc/panel_samples.txt          panel sample ids (column order)
  mhc/cohort_classII_gt.tsv      sample  n_sites  GT GT GT ... (tab-joined, phased a|b)
  relatives_ref/meta.json.gz     id -> breed label for panel samples
Outputs: mhc/panel_hap_groups.tsv, mhc/cohort_mhc.tsv, summary on stdout.
"""
import gzip, json, os, sys
import numpy as np

D = os.environ.get('D', '/u/project/pellegrini/gkislik/dogs')
# MHC_WIN=classII (default) or classI: file suffix + core block
WIN = os.environ.get('MHC_WIN', 'classII')
CORES = {'classII': (2330000, 2600000),     # DLA-DRB1 .. DLA-DQB1 (canFam4 chr12)
         'classI':  (1030000, 1175000)}     # DLA-88 / DLA-64
CORE = CORES[WIN]
CLUSTER_THR = 0.02             # haplotypes within 2% mismatches (core block) = same group
MIN_GROUP = 2                  # panel haplotypes needed to name a group


def load_sites():
    pos = [int(l.split('\t')[1]) for l in open(f'{D}/mhc/{WIN}_sites.tsv')]
    return np.array(pos)


def gt_to_haps(gt_fields):
    """list of 'a|b' strings -> (h1, h2) int8 arrays; missing/unphased -> -1"""
    n = len(gt_fields); h1 = np.full(n, -1, np.int8); h2 = np.full(n, -1, np.int8)
    for i, g in enumerate(gt_fields):
        if len(g) >= 3 and g[1] in '|/' and g[0] in '01' and g[2] in '01':
            h1[i] = int(g[0]); h2[i] = int(g[2])
    return h1, h2


def load_panel(pos):
    ids = [l.strip() for l in open(f'{D}/mhc/panel_samples.txt')]
    rows = [l.rstrip('\n').split('\t') for l in open(f'{D}/mhc/panel_{WIN}_gt.tsv')]
    assert len(rows) == len(pos), (len(rows), len(pos))
    gts = np.array([r[4:] for r in rows])            # sites x samples
    assert gts.shape[1] == len(ids)
    H = np.full((2 * len(ids), len(pos)), -1, np.int8)
    for j in range(len(ids)):
        h1, h2 = gt_to_haps(gts[:, j]); H[2 * j] = h1; H[2 * j + 1] = h2
    meta = json.load(gzip.open(f'{D}/relatives_ref/meta.json.gz', 'rt'))
    label = {s['id']: s['label'] for s in meta['samples']}
    breeds = [label.get(i, '?') for i in ids]
    return ids, breeds, H


def load_cohort(pos):
    ids, H = [], []
    for l in open(f'{D}/mhc/cohort_{WIN}_gt.tsv'):
        p = l.rstrip('\n').split('\t')
        if len(p) < 3: continue
        s, n, gts = p[0], int(p[1]), [g for g in p[2:] if g]
        if len(gts) != len(pos): print(f'  skip {s}: {len(gts)} GTs', file=sys.stderr); continue
        h1, h2 = gt_to_haps(gts); ids.append(s); H.append(h1); H.append(h2)
    return ids, np.array(H, np.int8)


def greedy_groups(Hc, thr):
    """cluster haplotypes (rows, core sites only) by Hamming fraction <= thr to the
    first member of each group (greedy, most-common-first). Returns group id per row."""
    n, m = Hc.shape
    keys = [Hc[i].tobytes() for i in range(n)]
    from collections import Counter
    cnt = Counter(keys)
    order = sorted(range(n), key=lambda i: (-cnt[keys[i]], keys[i]))
    centroids, cent_idx = [], []
    group = np.full(n, -1, int)
    for i in order:
        if group[i] >= 0: continue
        if centroids:
            C = np.array(centroids); d = (C != Hc[i]).mean(axis=1)
            j = int(d.argmin())
            if d[j] <= thr: group[i] = j; continue
        centroids.append(Hc[i].copy()); cent_idx.append(i); group[i] = len(centroids) - 1
        # absorb identical keys at once
    # second pass: assign everything to nearest centroid (fixes greedy order effects)
    C = np.array(centroids)
    for i in range(n):
        d = (C != Hc[i]).mean(axis=1); j = int(d.argmin())
        group[i] = j if d[j] <= thr else -1
    return group, C


def main():
    pos = load_sites()
    core = (pos >= CORE[0]) & (pos <= CORE[1])
    print(f'window sites {len(pos)}, core sites {core.sum()}')
    pids, pbreeds, PH = load_panel(pos)
    print(f'panel: {len(pids)} dogs, {len(set(pbreeds))} breeds')
    # panel haplotype groups on the core block
    pg, C = greedy_groups(PH[:, core], CLUSTER_THR)
    from collections import Counter, defaultdict
    gsize = Counter(pg)
    named = {g for g, c in gsize.items() if g >= 0 and c >= MIN_GROUP}
    print(f'panel haplotype groups: {len(set(pg)) - (1 if -1 in gsize else 0)} total, {len(named)} with >= {MIN_GROUP} haplotypes; '
          f'top 10 sizes {[c for g, c in gsize.most_common(10)]}')
    # per-group breed composition
    gb = defaultdict(Counter)
    for i, g in enumerate(pg): gb[g][pbreeds[i // 2]] += 1
    with open(f'{D}/mhc/panel_hap_groups_{WIN}.tsv', 'w') as f:
        f.write('group\tn_haplotypes\tn_breeds\ttop_breeds\n')
        for g, c in gsize.most_common():
            if g < 0: continue
            f.write(f'{g}\t{c}\t{len(gb[g])}\t' + '; '.join(f'{b} {k}' for b, k in gb[g].most_common(5)) + '\n')
    # panel per-dog stats: het rate, homozygous-group fraction, per-breed expectations
    pdiff = (PH[:, core][0::2] != PH[:, core][1::2]).mean(axis=1)
    phet = ((PH[0::2] != PH[1::2]) & (PH[0::2] >= 0) & (PH[1::2] >= 0)).mean(axis=1)
    phom = (pg[0::2] == pg[1::2]) & (pg[0::2] >= 0)
    breed_stats = defaultdict(list)
    for i, b in enumerate(pbreeds): breed_stats[b].append((phet[i], phom[i]))
    print(f'panel dogs: mean het {phet.mean():.3f}, homozygous class II haplotype {phom.mean()*100:.0f}%')
    worst = sorted(((np.mean([h for h, _ in v]), np.mean([m for _, m in v]), b, len(v)) for b, v in breed_stats.items() if len(v) >= 5))
    print('  least diverse breeds (mean het, %homozygous, n):', [(round(h, 3), int(m * 100), b, n) for h, m, b, n in worst[:8]])
    print('  most diverse breeds:', [(round(h, 3), int(m * 100), b, n) for h, m, b, n in worst[-6:]])

    if not os.path.exists(f'{D}/mhc/cohort_{WIN}_gt.tsv'):
        print('no cohort file yet'); return
    cids, CH = load_cohort(pos)
    print(f'cohort: {len(cids)} dogs')
    chet = ((CH[0::2] != CH[1::2]) & (CH[0::2] >= 0) & (CH[1::2] >= 0)).mean(axis=1)
    cdiff = (CH[:, core][0::2] != CH[:, core][1::2]).mean(axis=1)
    # assign cohort haplotypes to panel groups
    cg = np.full(CH.shape[0], -1, int)
    Cc = C
    for i in range(CH.shape[0]):
        d = (Cc != CH[i, core]).mean(axis=1); j = int(d.argmin())
        cg[i] = j if d[j] <= CLUSTER_THR else -1
    chom = (cg[0::2] == cg[1::2]) & (cg[0::2] >= 0)
    novel = (cg < 0).mean()
    print(f'cohort: mean het {chet.mean():.3f}, homozygous class II haplotype {chom.mean()*100:.0f}%, haplotypes not matching any panel group {novel*100:.1f}%')
    # depth / inbreeding / breed context
    depth, froh, top = {}, {}, {}
    for s in cids:
        r = f'{D}/results_prosper/{s}' if s.startswith('pk-') else f'{D}/results/{s.lower()}'
        try: depth[s] = json.load(open(f'{r}/qc_result.json'))['genome_mean_depth']
        except Exception: pass
        try: froh[s] = json.load(open(f'{r}/inbreeding_result.json'))['f_roh']
        except Exception: pass
        try:
            b = json.load(open(f'{r}/breed_result.json'))
            comp = b.get('breed_composition_global') or b['breed_composition']
            top[s] = (comp[0].get('breed_name') or comp[0]['breed'], float(comp[0]['proportion']))
        except Exception: pass
    order = np.argsort(chet); pct = np.empty(len(cids)); pct[order] = np.arange(len(cids)) / max(1, len(cids) - 1) * 100
    with open(f'{D}/mhc/cohort_mhc_{WIN}.tsv', 'w') as f:
        f.write('sample\thet_rate\thet_pct_cohort\tcore_hap_distance\tgroup1\tgroup2\thomozygous_classII\tdepth\tf_roh\ttop_breed\ttop_prop\n')
        for i, s in enumerate(cids):
            f.write(f'{s}\t{chet[i]:.4f}\t{pct[i]:.1f}\t{cdiff[i]:.4f}\t{cg[2*i]}\t{cg[2*i+1]}\t{int(chom[i])}\t{depth.get(s, "")}\t{froh.get(s, "")}\t{top.get(s, ("", ""))[0]}\t{top.get(s, ("", ""))[1]}\n')
    # bias checks
    ds = np.array([depth.get(s, np.nan) for s in cids]); fr = np.array([froh.get(s, np.nan) for s in cids])
    ok = ~np.isnan(ds)
    for lo, hi in [(0, 0.3), (0.3, 0.6), (0.6, 1), (1, 2), (2, 99)]:
        m = ok & (ds >= lo) & (ds < hi)
        if m.sum(): print(f'  depth {lo}-{hi}x: n={m.sum()} mean het {chet[m].mean():.3f} homozygous {chom[m].mean()*100:.0f}%')
    ok2 = ~np.isnan(fr)
    if ok2.sum() > 10: print(f'  corr(het_rate, F_ROH) = {np.corrcoef(chet[ok2], fr[ok2])[0,1]:.3f}; corr(core distance, F_ROH) = {np.corrcoef(cdiff[ok2], fr[ok2])[0,1]:.3f}')
    print('DONE')


if __name__ == '__main__':
    main()
