#!/usr/bin/env python3
"""
Stage 9e: MHC (DLA) diversity for one dog -> mhc_result.json.

  mhc_stage.py <imputed_dog10k.bcf> <mhc_reference.json> <out.json>
               [--breed breed_result.json] [--inbreeding inbreeding_result.json]

Per window (class II: DLA-DRA/DRB1/DQA1/DQB1; class I: DLA-88/DLA-64), from the
GLIMPSE2 phased genotypes at the Dog10K panel sites:
  het_rate               fraction of heterozygous sites in the window
  het_percentile         rank of het_rate among the processed cohort (reference)
  core_hap_distance      fraction of core-block sites at which the dog's two
                         phased haplotypes differ
  homozygous_haplotype   core_hap_distance <= hom_thr (0.02): the dog carries two
                         copies of essentially the same DLA haplotype
  breed_context          panel dogs of the dog's top breed: mean het, % homozygous
  vs_inbreeding          het expected from the dog's genome-wide F_ROH (cohort
                         regression) and the dog's residual in SD units
Validated 2026-09-13: imputed hets vs 5.5x reads at DP>=8 (Cosmo, class II):
sensitivity 1.00, precision 0.945; no depth bias 0.2-5x (analysis/mhc/).
"""
import json, sys
import numpy as np
import pysam


def haps(bcf, chrom, lo, hi):
    h1, h2 = [], []
    for rec in bcf.fetch(chrom, lo - 1, hi):
        gt = rec.samples[0].get('GT')
        if gt is None or len(gt) != 2 or gt[0] is None or gt[1] is None: continue
        h1.append(int(gt[0] > 0)); h2.append(int(gt[1] > 0))
    return np.array(h1, np.int8), np.array(h2, np.int8), [r.pos for r in bcf.fetch(chrom, lo - 1, hi)]


def main():
    args = sys.argv[1:]; bcf_path, ref_path, out_path = args[:3]
    def opt(n): return args[args.index(n) + 1] if n in args else None
    ref = json.load(open(ref_path)); hom_thr = ref['hom_thr']; chrom = ref['chrom']
    top_breed, top_prop, froh = None, None, None
    try:
        b = json.load(open(opt('--breed')))
        comp = b.get('breed_composition_global') or b['breed_composition']
        top_breed = comp[0].get('breed_name') or comp[0]['breed']; top_prop = float(comp[0]['proportion'])
    except Exception: pass
    try: froh = float(json.load(open(opt('--inbreeding')))['f_roh'])
    except Exception: pass

    bcf = pysam.VariantFile(bcf_path)
    res = {'chrom': chrom, 'hom_thr': hom_thr, 'windows': {}, 'top_breed': top_breed, 'top_breed_proportion': top_prop, 'f_roh': froh}
    for win, R in ref['windows'].items():
        (lo, hi), (clo, chi) = R['window'], R['core']
        h1, h2, pos = haps(bcf, chrom, lo, hi)
        pos = np.array(pos)
        if len(pos) < 0.5 * R['n_sites']:
            res['windows'][win] = {'status': 'insufficient_sites', 'n_sites': int(len(pos))}; continue
        het = float((h1 != h2).mean())
        core = (pos >= clo) & (pos <= chi)
        dist = float((h1[core] != h2[core]).mean()) if core.sum() else None
        sorted_het = np.array(R['cohort_het_sorted'])
        pct = float(np.searchsorted(sorted_het, het, side='right') / len(sorted_het) * 100)
        w = {'status': 'ok', 'n_sites': int(len(pos)), 'het_rate': round(het, 4), 'het_percentile': round(pct, 1),
             'core_hap_distance': round(dist, 4) if dist is not None else None,
             'homozygous_haplotype': bool(dist is not None and dist <= hom_thr),
             'cohort_mean_het': R['cohort_mean_het'], 'cohort_homozygous_frac': R['cohort_homozygous_frac']}
        if top_breed and top_breed in R['breeds'] and (top_prop or 0) >= 0.5:
            bc = R['breeds'][top_breed]
            w['breed_context'] = {'breed': top_breed, 'n_panel_dogs': bc['n'], 'mean_het': bc['mean_het'], 'homozygous_frac': bc['homozygous_frac']}
        if froh is not None:
            f = R['het_vs_froh']; exp = f['intercept'] + f['slope'] * froh
            w['vs_inbreeding'] = {'expected_het_from_froh': round(exp, 4), 'residual_sd': round((het - exp) / max(f['resid_sd'], 1e-6), 2)}
        res['windows'][win] = w

    c2 = res['windows'].get('classII', {})
    if c2.get('status') == 'ok':
        if c2['homozygous_haplotype']:
            level = 'low'
        elif c2['het_percentile'] < 25: level = 'below average'
        elif c2['het_percentile'] > 75: level = 'high'
        else: level = 'average'
        res['summary'] = {'level': level, 'classII_het_percentile': c2['het_percentile'], 'classII_homozygous': c2['homozygous_haplotype']}
    res['method'] = ('DLA (dog MHC) diversity from GLIMPSE2 phased genotypes at Dog10K panel sites on chr12 (class II: '
                     'DLA-DRA/DRB1/DQA1/DQB1 window; class I: DLA-88/DLA-64). het_rate = heterozygous fraction; a dog is '
                     'homozygous for the haplotype when its two phased core haplotypes differ at <= 2% of sites. Percentiles '
                     f'against {ref["windows"]["classII"]["cohort_n"]} ProsperK9 dogs; breed context from Dog10K panel dogs. '
                     'Imputed heterozygosity validated against 5.5x reads (sensitivity 1.00, precision 0.95).')
    json.dump(res, open(out_path, 'w'), indent=2)
    s = res.get('summary', {})
    print(f"mhc: classII het {c2.get('het_rate')} pct {c2.get('het_percentile')} homozygous {c2.get('homozygous_haplotype')} -> {s.get('level')}; "
          f"classI het {res['windows'].get('classI', {}).get('het_rate')}")


if __name__ == '__main__':
    main()
