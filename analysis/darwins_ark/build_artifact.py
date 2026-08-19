#!/usr/bin/env python3
"""Assemble the Darwin's Ark individual-trait artifact for Stage 11.

Outputs into ../reference_json/darwins_ark/:
  wts_<trait>.tsv.gz    chr pos ea oa beta af  (platform-filtered, p<=0.1)
  manifest.json.gz      per-trait reference distribution (sorted cohort PRS on
                        the SAME sites), truth mean/sd, category anchors,
                        held-out cv_r, h2 (Morrill 2022), ancestry-blend
                        coefficients + per-breed trait means keyed by OUR
                        230-breed panel labels
"""
import gzip, json, os, shutil, numpy as np

OUT='../reference_json/darwins_ark'
os.makedirs(OUT, exist_ok=True)
FACTORS=['human_sociability','arousal_level','toy_motor_patterns','biddability',
         'agonistic_threshold','dog_sociability','environmental_engagement','proximity_seeking']
PHYS=['white_fur','curly_tail','ear_shape','fur_length','fur_texture']
DISPLAY={'human_sociability':'Human Sociability','arousal_level':'Arousal Level',
 'toy_motor_patterns':'Toy-directed Motor Patterns','biddability':'Biddability',
 'agonistic_threshold':'Agonistic Threshold','dog_sociability':'Dog Sociability',
 'environmental_engagement':'Environmental Engagement','proximity_seeking':'Proximity Seeking',
 'size':'Body Size','white_fur':'White Fur Amount','curly_tail':'Curly Tail',
 'ear_shape':'Ear Shape','fur_length':'Fur Length','fur_texture':'Fur Texture'}
H2={'human_sociability':0.41,'arousal_level':0.23,'toy_motor_patterns':0.30,'biddability':0.27,
    'agonistic_threshold':0.10,'dog_sociability':0.15,'environmental_engagement':0.18,
    'proximity_seeking':0.17,'size':0.93,'white_fur':0.46,'curly_tail':0.90,'ear_shape':0.99,
    'fur_length':0.99,'fur_texture':0.39}
ANCHORS={'curly_tail':[[-0.601,'not curly'],[1.6635,'curly']],
 'fur_length':[[-0.881,'short'],[0.4755,'medium'],[1.8321,'long']],
 'fur_texture':[[-0.426,'soft'],[2.345,'rough/wiry']],
 'ear_shape':[[-1.4,'floppy'],[1.35,'erect']],
 'white_fur':[[-1.5,'lots of white'],[0.5,'little white']]}
CVR={'size':0.774,'human_sociability':0.212,'biddability':0.241,'toy_motor_patterns':0.222,
     'arousal_level':0.156,'dog_sociability':0.086,'environmental_engagement':0.096,
     'proximity_seeking':0.138,'agonistic_threshold':0.080}

blendinfo=json.load(open('da_blend_ancestry.json'))

# truth stats
cols=[l.rstrip('\n').split('\t') for l in open('gemma/columns.tsv')]
pheno=[l.rstrip('\n').split('\t') for l in open('gemma/pheno.tsv')]
pheno_phys=[l.rstrip('\n').split('\t') for l in open('gemma/pheno_phys.tsv')]
phys_order=['white_fur','curly_tail','ear_shape','fur_length','fur_texture']
def truth_stats(trait):
    if trait in phys_order:
        i=phys_order.index(trait)
        v=[float(r[i]) for r in pheno_phys if r[i] not in ('NA','')]
    else:
        col=next(int(c) for t,f,c in cols if t==trait and f=='full')
        v=[float(r[col-1]) for r in pheno if r[col-1]!='NA']
    return float(np.mean(v)), float(np.std(v)), len(v)

manifest={'meta':{'source':"Darwin's Ark (Morrill et al. 2022, Science abk0639; Dryad CC0)",
                  'model':'GEMMA -lmm 1, GRM + sex/data-type(+age,size) covariates, 2,155 dogs',
                  'policy':'p<=0.1 weights, platform-filtered to Dog10K-imputable sites',
                  'reference':'cohort PRS distribution on the same filtered sites',
                  'built':'2026-08-19'},
          'traits':{}}
for trait in FACTORS+['size']+PHYS:
    src=f'artifact_wts_{trait}.tsv.gz'
    shutil.copy(src, f'{OUT}/wts_{trait}.tsv.gz')
    ref=np.load(f'artifact_ref_{trait}.npy')
    mu,sd,n=truth_stats(trait)
    t={'display':DISPLAY[trait],
       'kind':'factor' if trait in FACTORS else ('physical' if trait in PHYS else 'size'),
       'n_sites':sum(1 for _ in gzip.open(f'{OUT}/wts_{trait}.tsv.gz','rt')),
       'h2_snp':H2[trait],
       'truth_mean':round(mu,4),'truth_sd':round(sd,4),'n_pheno':n,
       'ref_prs_sorted':[round(float(x),3) for x in ref]}
    if trait in CVR: t['cv_r_heldout']=CVR[trait]
    if trait in ANCHORS: t['anchors']=ANCHORS[trait]
    if trait in blendinfo['blend']: t['ancestry_blend']=blendinfo['blend'][trait]
    manifest['traits'][trait]=t
manifest['breed_trait_means']=blendinfo['breed_trait_means']
with gzip.open(f'{OUT}/manifest.json.gz','wt') as fh:
    json.dump(manifest,fh,separators=(',',':'))
tot=sum(os.path.getsize(f'{OUT}/{f}') for f in os.listdir(OUT))
print(f'artifact: {len(manifest["traits"])} traits, {tot/1e6:.1f} MB in {OUT}')
