"""Run the ASIP caller on purebred cohort dogs (top breed >= 85%) from breeds
whose A-locus pattern is fixed or near-fixed, for comparison with breed standards
and Bannasch et al. 2021 Supplementary Table 6."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from asip_locus import call_asip
EXPECT = {
    'PUG': 'ay', 'FRENCH_BULLDOG': 'ay', 'BOXER': 'ay', 'SHIBA_INU': 'ay', 'BULLMASTIFF': 'ay',
    'PEMBROKE_WELSH_CORGI': 'ay', 'CHIHUAHUA': 'ay', 'BASENJI': 'ay',
    'ROTTWEILER': 'at', 'DOBERMAN_PINSCHER': 'at', 'BEAGLE': 'at', 'BERNESE_MOUNTAIN': 'at',
    'YORKSHIRE_TERRIER': 'at', 'GERMAN_SHEPHERD': 'at', 'MINIATURE_PINSCHER': 'at', 'AIREDALE_TERRIER': 'at',
    'SIBERIAN_HUSKY': 'aw', 'ALASKAN_MALAMUTE': 'aw', 'KEESHOND': 'aw', 'NORWEGIAN_ELKHOUND': 'aw',
}
H, LIST, PER = sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 8
n = {}
for line in open(LIST):
    pk, bc, pub, work = line.rstrip('\n').split('\t')[:4]
    try:
        b = json.load(open(f'{pub}/breed_result.json'))['breed_composition'][0]
        q = json.load(open(f'{pub}/qc_result.json'))
    except Exception: continue
    code = b['breed']
    if code not in EXPECT or b['proportion'] < 0.85 or n.get(code, 0) >= PER: continue
    bam = f'{work}/markdup.bam'
    if not os.path.exists(bam): continue
    n[code] = n.get(code, 0) + 1
    r = call_asip(bam, H)
    print('\t'.join(map(str, [pk, code, EXPECT[code], round(b['proportion'], 2), q.get('genome_mean_depth'),
          '/'.join(r['alleles']), r['confidence'], json.dumps(r['hcp_reads']), json.dumps(r['tag_reads']),
          json.dumps(r['r96c_reads']), json.dumps(r['vp_reads']), list(r['class_post'].items())[0]])), flush=True)
