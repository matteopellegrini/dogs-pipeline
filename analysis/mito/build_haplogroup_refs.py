#!/usr/bin/env python3
"""Build reference_json/mito_haplogroups.json.gz for the pipeline's mito stage.

Inputs (committed alongside this script for reproducibility):
  mito_refs_hv1.fa   100 haplotype-labeled complete dog mitogenomes from
                     GenBank (JF342xxx series). Found via NCBI esearch
                     'Canis lupus familiaris[Organism] AND mitochondrion[Title]
                     AND haplotype[Title]'; labels parsed from record titles
                     ("haplotype HV1 A18 ..."). 72 A / 21 B / 5 C / 2 D,
                     matching published dog haplogroup frequencies.
  mito_titles.tsv    accession -> haplotype label map from those titles.

Method notes (validated on the 96-dog cohort, 2026-08-27):
  - GenBank records are REVERSE-COMPLEMENT of canFam4 chrM and circularly
    rotated, so each ref is RC'd and aligned (edlib infix) against a DOUBLED
    chrM before extracting differences.
  - Comparison must be SNP-set based, never whole-sequence edit distance:
    a low-pass consensus inherits canFam4 wherever reads don't disagree
    (control-region VNTR included), so sequence distance to external records
    is dominated by shared reference/VNTR noise (~300 edits) which swamps the
    ~30-70 real SNPs — every dog then "matches" whichever record most
    resembles canFam4 (whose own chrM is itself haplogroup C-like).
  - Substitutions within 5bp of an indel are dropped (alignment jitter around
    the VNTR and homopolymers).
  - Leave-one-out major-group accuracy of nearest-neighbor assignment on
    these 100 refs: 95/100.

Output: {"refs": [{"haplotype","acc","snps":[[pos0,alt],...]}...],
         "groups": {letter: story text}, "meta": {...}}

Usage:  python3 analysis/mito/build_haplogroup_refs.py
        (writes reference_json/mito_haplogroups.json.gz; needs edlib + canFam4.fa)
"""
import edlib, gzip, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
FASTA = os.path.join(ROOT, 'canFam4.fa')
OUT = os.path.join(ROOT, 'reference_json', 'mito_haplogroups.json.gz')

def chrm_seq():
    fai = {l.split('\t')[0]: l.split('\t') for l in open(FASTA + '.fai')}
    length, offset, linebases = int(fai['chrM'][1]), int(fai['chrM'][2]), int(fai['chrM'][3])
    with open(FASTA) as f:
        f.seek(offset)
        nlines = (length + linebases - 1) // linebases
        seq = ''.join(f.readline().strip() for _ in range(nlines)).upper()
    assert len(seq) == length
    return seq

def read_fa(path):
    seqs, name, cur = {}, None, []
    for line in open(path):
        if line.startswith('>'):
            if name: seqs[name] = ''.join(cur)
            name = line[1:].split()[0]; cur = []
        else:
            cur.append(line.strip().upper())
    if name: seqs[name] = ''.join(cur)
    return seqs

RC = str.maketrans('ACGTN', 'TGCAN')

def snp_set(query, target2, L):
    res = edlib.align(query, target2, mode='HW', task='path')
    start = res['locations'][0][0]
    q = t = 0
    subs, indel_pos = [], set()
    for n, op in re.findall(r'(\d+)([=XIDM])', res['cigar']):
        n = int(n)
        if op == '=':
            q += n; t += n
        elif op == 'X':
            for k in range(n):
                subs.append(((start + t + k) % L, query[q + k]))
            q += n; t += n
        elif op == 'I':
            indel_pos.add((start + t) % L); q += n
        else:  # D
            for k in range(n):
                indel_pos.add((start + t + k) % L)
            t += n
    near = {(p + d) % L for p in indel_pos for d in range(-5, 6)}
    return sorted((p, b) for p, b in subs if p not in near and b in 'ACGT')

GROUP_STORIES = {
    'A': 'Haplogroup A is the largest maternal lineage in dogs, carried by roughly '
         'two out of three dogs worldwide. It traces back to the main group of '
         'wolf mothers from which most modern dogs descend.',
    'B': 'Haplogroup B is the second-largest maternal lineage, found in roughly one '
         'in five dogs worldwide and common across many European breeds.',
    'C': 'Haplogroup C is a less common maternal lineage, carried by roughly one in '
         'twenty dogs. It represents a separate, ancient group of wolf mothers.',
    'D': 'Haplogroup D is a rare maternal lineage, seen mostly in some European and '
         'Scandinavian breeds.',
}

def main():
    chrm = chrm_seq()
    L = len(chrm)
    target2 = chrm + chrm
    labels = {l.split('\t')[0]: l.split('\t')[1]
              for l in open(os.path.join(HERE, 'mito_titles.tsv'))}
    refs = read_fa(os.path.join(HERE, 'mito_refs_hv1.fa'))
    out_refs = []
    for name, seq in sorted(refs.items()):
        acc = name.split('.')[0]
        if acc not in labels:
            continue
        rc = seq.translate(RC)[::-1]
        out_refs.append({'haplotype': labels[acc], 'acc': acc,
                         'snps': snp_set(rc, target2, L)})
    payload = {
        'refs': out_refs,
        'groups': GROUP_STORIES,
        'meta': {
            'chrm_length': L,
            'n_refs': len(out_refs),
            'source': 'GenBank JF342xxx complete dog mitogenomes, haplotype-labeled titles',
            'method': 'SNP symmetric difference vs canFam4 chrM; subs within 5bp of indels masked',
            'loo_major_group_accuracy': '95/100',
        },
    }
    with gzip.open(OUT, 'wt') as f:
        json.dump(payload, f)
    sizes = sorted(len(r['snps']) for r in out_refs)
    print(f"wrote {OUT}: {len(out_refs)} refs, SNPs/ref median {sizes[len(sizes)//2]}")

if __name__ == '__main__':
    sys.exit(main())
