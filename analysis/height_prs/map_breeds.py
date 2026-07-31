"""Map Parker breed codes -> AKC breed names in the height dataset."""
import csv, re, ast, sys, json

PIPE = "/Users/matteopellegrini/Downloads/dogs/run_dog_pipeline.sh"
FAM  = "/Users/matteopellegrini/Downloads/dogs/COSMO/analysis/cosmo_parker_full.fam"
AKC  = "/private/tmp/claude-501/-Users-matteopellegrini-Downloads-dogs/0fb26461-2afd-4142-936c-a58295d35959/scratchpad/height/akc.csv"

# ── pull PARKER_TO_AKC out of the pipeline so we don't re-type 150 entries ──
src = open(PIPE).read()
start = src.index("PARKER_TO_AKC = {")
body = src[start + len("PARKER_TO_AKC = "):]
depth = 0
for i, ch in enumerate(body):
    if ch == '{': depth += 1
    elif ch == '}':
        depth -= 1
        if depth == 0:
            body = body[:i+1]; break
body = re.sub(r'#.*', '', body)
PARKER_TO_AKC = ast.literal_eval(body)

akc_rows = list(csv.DictReader(open(AKC)))
heights = {}
for r in akc_rows:
    name = r[''].strip()
    try:
        lo, hi = float(r['min_height']), float(r['max_height'])
    except (ValueError, TypeError):
        continue
    heights[name] = (lo + hi) / 2.0     # cm
akc_names = set(heights)

# kkakey plural/group style -> AKC singular style
def kkakey_to_akc(n):
    m = re.match(r'^(.+?)s \((.+)\)$', n)           # "Retrievers (Golden)" -> "Golden Retriever"
    if m: return f"{m.group(2)} {m.group(1)}"
    irregular = {
        'Keeshonden':'Keeshond', 'Pulik':'Puli', 'Pumik':'Pumi',
        'Komondorok':'Komondor', 'Kuvaszok':'Kuvasz', 'Maltese':'Maltese',
        'Pekingese':'Pekingese', 'Havanese':'Havanese', 'Shih Tzu':'Shih Tzu',
        'Chinese Shar-Pei':'Chinese Shar-Pei', 'Cane Corso':'Cane Corso',
        'Belgian Tervuren':'Belgian Tervuren', 'Shiba Inu':'Shiba Inu',
        'Chinese Crested':'Chinese Crested', 'Xoloitzcuintli':'Xoloitzcuintli',
        'Coton de Tulear':'Coton de Tulear', 'Great Pyrenees':'Great Pyrenees',
        'Japanese Chin':'Japanese Chin', 'Finnish Spitz':'Finnish Spitz',
        'Petits Bassets Griffons Vendeens':'Petit Basset Griffon Vendeen',
        'Dogues de Bordeaux':'Dogue de Bordeaux',
        'Bouviers des Flandres':'Bouvier des Flandres',
        'Bichons Frises':'Bichon Frise',
        'Spinoni Italiani':'Spinone Italiano',
        'Wirehaired Pointing Griffons':'Wirehaired Pointing Griffon',
        'Russell Terriers':'Russell Terrier',
        'Belgian Malinois':'Belgian Malinois',
    }
    if n in irregular: return irregular[n]
    if n.endswith('ies'):  return n[:-3] + 'y'
    if n.endswith('es') and n[:-2] in akc_names: return n[:-2]
    if n.endswith('s'):    return n[:-1]
    return n

# Size varieties the pipeline collapses but height cannot
OVERRIDE = {
    'SPOO': 'Poodle (Standard)', 'MPOO': 'Poodle (Miniature)', 'TPOO': 'Poodle (Toy)',
    'ACKR': 'Cocker Spaniel',            # AKC "Cocker Spaniel" is the American one
    'COOK': 'Cocker Spaniel',
    'BORD': 'Border Collie', 'COLL': 'Collie',
    'PBGV': 'Petit Basset Griffon Vendéen',
    'TYFX': 'Smooth Fox Terrier', 'WFOX': 'Wire Fox Terrier',
    # Real AKC breeds the pipeline's behavioural-trait map happened to omit
    'STBD': 'Saint Bernard', 'NELK': 'Norwegian Elkhound', 'NOWT': 'Norwich Terrier',
    'REDB': 'Redbone Coonhound', 'EURA': 'Eurasier', 'CIRN_Italy': 'Cirneco dell’Etna',
    'KELP': 'Australian Kelpie', 'INCA': 'Peruvian Inca Orchid',
    'FOXH': 'American Foxhound',   # Parker code ambiguous vs English; heights differ ~2cm
}

# Excluded on purpose:
#   MNTY  Manchester Terrier — Standard vs Toy differ ~10cm, code doesn't say which
#   LMUN/GDJK  unresolved codes
#   *_Italy / XIGO_China  landrace & village dogs with no AKC standard

# Not AKC breeds — wolves, village/landrace dogs, and the query sample itself.
EXCLUDE_PREFIX = ('WOLF',)
EXCLUDE = {'COSMO'}

code_to_breed, unmatched = {}, []
codes = sorted({l.split()[0] for l in open(FAM)})
for code in codes:
    if code in EXCLUDE or code.startswith(EXCLUDE_PREFIX):
        continue
    if code in OVERRIDE:
        code_to_breed[code] = OVERRIDE[code]; continue
    kk = PARKER_TO_AKC.get(code)
    if not kk:
        unmatched.append((code, None)); continue
    cand = kkakey_to_akc(kk)
    if cand in akc_names:
        code_to_breed[code] = cand
    else:
        hit = next((n for n in akc_names if n.lower() == cand.lower()), None)
        if hit: code_to_breed[code] = hit
        else:   unmatched.append((code, f"{kk} -> {cand}"))

print(f"Parker codes: {len(codes)}   mapped: {len(code_to_breed)}   unmatched: {len(unmatched)}")
print("\nUnmatched:")
for c, why in unmatched: print(f"  {c:<14} {why}")

json.dump({'code_to_breed': code_to_breed,
           'height_cm': {b: heights[b] for b in set(code_to_breed.values())}},
          open('/private/tmp/claude-501/-Users-matteopellegrini-Downloads-dogs/0fb26461-2afd-4142-936c-a58295d35959/scratchpad/height/mapping.json','w'), indent=1)
