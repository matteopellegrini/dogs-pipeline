"""Dry run of the A-locus caller over a rerun list (read-only): one TSV row per
kit with the new call and the A call currently published.
usage: cohort_scan.py asip_haplotypes.json list.tsv start end > out.tsv"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from asip_locus import call_asip
H, LIST, a, b = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
for i, line in enumerate(open(LIST), 1):
    if not a <= i <= b: continue
    pk, bc, pub, work = line.rstrip('\n').split('\t')[:4]
    bam = f'{work}/markdup.bam'
    if not os.path.exists(bam): print(f'{pk}\tNO_BAM'); continue
    try:
        old = json.load(open(f'{pub}/coat_color.json'))['loci']['A']['predicted_alleles']
        dp = json.load(open(f'{pub}/qc_result.json')).get('genome_mean_depth')
        top = json.load(open(f'{pub}/breed_result.json'))['breed_composition'][0]
    except Exception as e:
        old, dp, top = None, None, {}
    r = call_asip(bam, H)
    print('\t'.join(map(str, [pk, dp, top.get('breed'), round(top.get('proportion', 0), 2), '/'.join(old or []),
                              '/'.join(r['alleles']), r['confidence'], r['informative_reads'], r['ay_support'],
                              json.dumps(r['r96c_reads']), json.dumps(r['vp_reads'])])), flush=True)
