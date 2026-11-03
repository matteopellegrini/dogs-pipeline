#!/usr/bin/env python3
"""Panel v3, step 2a: which frozen small variants have carriers in our cohort, and in which breeds.

Input is the cohort matrix built by cluster/panel_v3_merge.sh from the per-kit read genotypes
(work_panel_v3/cohort_reads.tsv.gz: CHROM POS REF ALT then one AD field per kit) plus the kit manifest
(pk, barcode, pub dir, work dir) for each kit's top breed.

Per site, with a per-read error rate EPS:
  * cohort allele frequency by EM over all kits' read counts (diploid, HWE, binomial reads);
  * per-kit carrier posterior P(>=1 alt copy | reads, cohort AF);
  * kits with alt reads, total alt reads, expected carriers (sum of posteriors), carrier breeds.
A site whose expected carrier count is ~0 cannot be imputed from the cohort; it stays read-only.

  python3 panel_v3_carriers.py cohort_reads.tsv.gz manifest.tsv sites.tsv OUT_PREFIX
"""
import gzip, json, os, sys, collections
import numpy as np

EPS = 0.005  # per-read error after Q20/MAPQ20 filtering

def load_matrix(path):
    with gzip.open(path, 'rt') as f:
        header = f.readline().rstrip('\n').split('\t')
        samples = header[4:]
        sites, ref_c, alt_c = [], [], []
        for line in f:
            p = line.rstrip('\n').split('\t')
            sites.append(tuple(p[:4]))
            r, a = [], []
            for ad in p[4:]:
                # AD is "." (kit absent), "7" (hom-ref record, no ALT) or "ref,alt[,alt2]" ("." for unseen alleles)
                v = [0 if x == '.' else int(x) for x in ad.split(',')]
                r.append(v[0]); a.append(sum(v[1:]))
            ref_c.append(r); alt_c.append(a)
    return samples, sites, np.array(ref_c, dtype=np.int32), np.array(alt_c, dtype=np.int32)

def genotype_lik(r, a):
    """P(reads | g) for g = 0, 1, 2 alt copies; arrays of shape (n,) -> (n, 3)."""
    p_alt = np.array([EPS, 0.5, 1 - EPS])
    # log binomial without the constant (identical across genotypes)
    ll = a[:, None] * np.log(p_alt)[None, :] + r[:, None] * np.log(1 - p_alt)[None, :]
    ll -= ll.max(axis=1, keepdims=True)
    return np.exp(ll)

def em_af(lik, iters=200):
    q = 0.01
    for _ in range(iters):
        prior = np.array([(1 - q) ** 2, 2 * q * (1 - q), q ** 2])
        post = lik * prior
        post /= post.sum(axis=1, keepdims=True)
        q_new = (post[:, 1] + 2 * post[:, 2]).sum() / (2 * len(post))
        if abs(q_new - q) < 1e-9:
            break
        q = q_new
    return q, post

def top_breed(pub):
    try:
        d = json.load(open(os.path.join(pub, 'breed_result.json')))
        b = d['breed_composition'][0]
        return f"{b.get('breed_name', b['breed'])} {round(100 * b['proportion'])}%"
    except Exception:
        return 'unknown'

def main(matrix, manifest, sites_tsv, out):
    samples, sites, R, A = load_matrix(matrix)
    pubs = {}
    for line in open(manifest):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3:
            pubs[p[0]] = p[2]
    breeds = [top_breed(pubs.get(s, '')) for s in samples]
    meta = {}
    for line in open(sites_tsv):
        p = line.rstrip('\n').split('\t')
        if p[0] == 'site_id':
            cols = p; continue
        row = dict(zip(cols, p))
        meta[(row['chrom'], row['pos'])] = row
    with open(out + '.sites.tsv', 'w') as fo, open(out + '.carriers.tsv', 'w') as fc:
        fo.write('chrom\tpos\tref\talt\tgene\tcondition\treport_class\tkits_covered\tkits_alt_reads\talt_reads'
                 '\ttotal_reads\tcohort_af\texpected_carriers\tconfident_carriers\tcarrier_breeds\n')
        fc.write('chrom\tpos\tgene\tkit\tbreed\tref_reads\talt_reads\tp_carrier\n')
        for i, (c, pos, ref, alt) in enumerate(sites):
            r, a = R[i], A[i]
            q, post = em_af(genotype_lik(r, a))
            pc = post[:, 1] + post[:, 2]
            m = meta.get((c, pos), {})
            idx = np.where((a > 0) | (pc > 0.5))[0]
            bc = collections.Counter(breeds[j].rsplit(' ', 1)[0] for j in idx if pc[j] > 0.5)
            fo.write('\t'.join(map(str, [c, pos, ref, alt, m.get('gene', ''), m.get('condition', ''),
                     m.get('report_class', ''), int(((r + a) > 0).sum()), int((a > 0).sum()), int(a.sum()),
                     int((r + a).sum()), f'{q:.5f}', f'{pc.sum():.1f}', int((pc > 0.9).sum()),
                     '; '.join(f'{k} x{v}' for k, v in bc.most_common(6))])) + '\n')
            for j in idx:
                fc.write('\t'.join(map(str, [c, pos, m.get('gene', ''), samples[j], breeds[j],
                         int(r[j]), int(a[j]), f'{pc[j]:.3f}'])) + '\n')
    print(f'{len(sites)} sites x {len(samples)} kits -> {out}.sites.tsv, {out}.carriers.tsv')

if __name__ == '__main__':
    main(*sys.argv[1:5])
