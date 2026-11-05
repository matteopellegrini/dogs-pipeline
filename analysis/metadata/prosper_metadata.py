#!/usr/bin/env python3
"""One metadata row per ProsperKits sample, keyed by the original FASTQ filename.

Joins, per pk-<kit> sample in the batch-1 and batch-2 sample sheets:
  * the archive FASTQ files the sample was built from (targets of the fastq_prosper symlinks);
  * customer-entered age and weight from prosperKitAgeInfo (free text, matched by kit number;
    entered at order time in 2024, so age is age at swab), plus parsed age_years / weight_lbs;
  * mean depth and QC status (qc_result.json), DNA sex (cnv50_result.json) and our breed estimate
    (breed_result.json, the published breed panel v2.3 composition).
DT-prefixed files are targeted bisulfite libraries, not WGS: they are left out of the file list, and samples
built only from them are skipped.

  python3 prosper_metadata.py AGE_INFO OUT.tsv SHEET [SHEET ...]
"""
import json, os, re, sys

HEDGE = re.compile(r'^(~|\+/-|about|approx\.?|approximately|around)\s*')

def parse_age(s):
    s = HEDGE.sub('', s.strip().lower()).replace('ish', '')
    m = re.match(r'^(\d+(?:\.\d+)?)(?![\d.])\s*(years?|yrs?|yo|y|months?|mon\w*|mths?|mos?|mo|m|weeks?|wks?|w)?(?![a-z])', s)
    if not m:
        return ''
    v, unit = float(m.group(1)), (m.group(2) or 'y')
    v = v / 12 if unit.startswith('m') else v / 52 if unit.startswith('w') else v
    return f'{v:.2f}' if v < 30 else ''

def parse_weight(s):
    s = HEDGE.sub('', s.strip().lower())
    m = re.match(r'^(\d+(?:\.\d+)?)(?:\s*-\s*(\d+(?:\.\d+)?))?', s)
    if not m:
        return ''
    v = float(m.group(1)) if not m.group(2) else (float(m.group(1)) + float(m.group(2))) / 2
    if 'kg' in s:
        v *= 2.20462
    return f'{v:.1f}' if 0 < v < 300 else ''

BISULFITE = re.compile(r'(^|--)DT\d')

def main(age_info, out, *sheets):
    info = {}
    for line in open(age_info):
        p = line.rstrip('\n').split('\t')
        if len(p) < 4 or not p[0].isdigit():
            continue
        kit, age, weight = p[0], p[1].strip(), p[2].strip()
        cur = info.setdefault(kit, ['', ''])
        if age and age.upper() != 'NO AGE' and not cur[0]:
            cur[0] = age
        if weight and not cur[1]:
            cur[1] = weight
    cols = ['fastq_r1', 'fastq_files', 'kit', 'sample_id', 'age_reported', 'age_years', 'weight_reported',
            'weight_lbs', 'sex_dna', 'mean_depth', 'qc_status', 'top_breed', 'top_breed_pct', 'breed_composition']
    n = skipped = 0
    seen = set()
    with open(out, 'w') as fo:
        fo.write('\t'.join(cols) + '\n')
        for sheet in sheets:
            for line in open(sheet):
                p = line.rstrip('\n').split('\t')
                if p[0] == 'sample_id' or not p[0].startswith('pk-') or p[0] in seen:
                    continue
                sid, fq, pub = p[0], p[1], p[5]
                seen.add(sid)
                files = sorted(os.path.basename(os.path.realpath(os.path.join(fq, f)))
                               for f in os.listdir(fq)) if os.path.isdir(fq) else []
                wgs = [f for f in files if not BISULFITE.search(f)]
                if files and not wgs:
                    skipped += 1
                    continue
                files = wgs
                r1 = [f for f in files if re.search(r'_R1(_\d+)?\.f(ast)?q\.gz$', f)]
                # primary = the plain library; -TR / -BR top-ups sort after it
                r1.sort(key=lambda f: (bool(re.search(r'-(TR|BR)_', f)), f))
                kit = sid[3:]
                age, weight = info.get(kit, ['', ''])
                row = dict(fastq_r1=r1[0] if r1 else '', fastq_files=';'.join(files), kit=kit, sample_id=sid,
                           age_reported=age, age_years=parse_age(age), weight_reported=weight,
                           weight_lbs=parse_weight(weight))
                try:
                    q = json.load(open(os.path.join(pub, 'qc_result.json')))
                    row.update(mean_depth=q.get('genome_mean_depth', ''), qc_status=q.get('qc_status', ''))
                except Exception:
                    row.update(qc_status='NO_RESULTS')
                try:
                    row['sex_dna'] = json.load(open(os.path.join(pub, 'cnv50_result.json'))).get('sex_inferred', '')
                except Exception:
                    pass
                try:
                    bc = json.load(open(os.path.join(pub, 'breed_result.json')))['breed_composition']
                    bc = [b for b in bc if b['proportion'] >= 0.005]
                    row.update(top_breed=bc[0].get('breed_name', bc[0]['breed']),
                               top_breed_pct=round(100 * bc[0]['proportion'], 1),
                               breed_composition='; '.join(f"{b.get('breed_name', b['breed'])} {100 * b['proportion']:.1f}%"
                                                           for b in bc))
                except Exception:
                    pass
                fo.write('\t'.join(str(row.get(c, '')) for c in cols) + '\n')
                n += 1
    print(f'{n} samples -> {out} ({skipped} bisulfite-only samples skipped)')

if __name__ == '__main__':
    main(*sys.argv[1:])
