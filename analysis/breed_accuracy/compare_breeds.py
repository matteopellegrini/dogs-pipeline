#!/usr/bin/env python3
"""
Compare SCOPE breed predictions against owner-reported breeds.

    python3 compare_breeds.py reported.tsv predictions.tsv parker_names.json

  reported.tsv     sample_no, name, reported_breed, sex
  predictions.tsv  sample, top1, p1, top2, p2, top3, p3   (breed CODES)
  parker_names.json  {code: full name}

Owner reports are free text, so most of this file is about matching them to the
Parker vocabulary honestly. Three things make a naive comparison misleading:

  Misspellings and abbreviations. "Australian Shepard", "Pomerarian",
  "St Bernard" are all correct reports of breeds the model got right. Fuzzy
  matching against the Parker vocabulary handles them; exact matching scores
  them as model errors.

  Size prefixes are breed-defining, not descriptive. Parker lists Miniature,
  Standard and Giant Schnauzer separately, so stripping "Miniature" as a
  descriptor turns a correct call into a miss. Colours ("chocolate", "black")
  genuinely are descriptive and are stripped.

  "Mix" means the owner is telling you it is NOT purebred. Scoring
  "Beagle Mix" as a failed Beagle prediction measures the wrong thing. These
  are reported separately, where the question is whether the named breed shows
  up at all, not whether it wins.

Designer crosses (Labradoodle, Cavachon) have no Parker label by construction.
They are scored on whether BOTH parent breeds appear, which is the only
prediction that could be right.
"""
import difflib
import json
import re
import sys

# Colours and marketing words carry no breed information. Size prefixes are
# deliberately NOT here: Parker treats them as distinct breeds.
DESCRIPTORS = {'CHOCOLATE', 'BLACK', 'YELLOW', 'WHITE', 'CREAM', 'RED', 'BLUE',
               'MERLE', 'GOLDEN', 'PURE', 'BRED', 'PUREBRED', 'DOG', 'PUPPY',
               'LILLAC', 'LILAC', 'ISABELLE'}

# Colloquial names owners use that are not Parker labels.
SYNONYMS = {
    'LAB': 'LABRADOR RETRIEVER', 'DOXIE': 'DACHSHUND', 'GSD': 'GERMAN SHEPHERD DOG',
    'ST': 'SAINT', 'SHEPARD': 'SHEPHERD', 'PIT BALL': 'PIT BULL',
    'PIT BULL': 'AMERICAN STAFFORDSHIRE TERRIER', 'PITBULL': 'AMERICAN STAFFORDSHIRE TERRIER',
}
# Designer crosses -> the two parent breeds the model should actually report.
CROSSES = {
    'LABRADOODLE': ('Labrador Retriever', 'Poodle'),
    'GOLDENDOODLE': ('Golden Retriever', 'Poodle'),
    'GOLDEN DOODLE': ('Golden Retriever', 'Poodle'),
    'CAVACHON': ('Cavalier King Charles Spaniel', 'Bichon Frise'),
    'COCKAPOO': ('American Cocker Spaniel', 'Poodle'),
    'PUGGLE': ('Pug', 'Beagle'),
}
MIX_RE = re.compile(r'\bMIX(ED)?\b', re.I)
SPLIT_RE = re.compile(r'[/+,&]| AND | X |-(?=[A-Za-z])', re.I)
MATCH_MIN = 0.82         # fuzzy ratio needed when tokens do not nest


def norm(s):
    s = re.sub(r'[^A-Za-z]+', ' ', s or '').upper()
    return re.sub(r'\s+', ' ', s).strip()


def tokens(s):
    t = norm(s)
    for k, v in SYNONYMS.items():
        t = re.sub(rf'\b{k}\b', v, t)
    return [w for w in t.split() if w not in DESCRIPTORS]


def canonical(text, vocab):
    """Best Parker breed for a free-text fragment -> (name, score)."""
    t = ' '.join(tokens(text))
    if not t:
        return None, 0.0
    best, score = None, 0.0
    for c in vocab:
        ct = ' '.join(tokens(c))
        if not ct:
            continue
        a, b = set(t.split()), set(ct.split())
        r = difflib.SequenceMatcher(None, t, ct).ratio()
        if not (a & b):
            # Requiring a shared word stops "Caucasian Shepherd" matching
            # "Australian Shepherd". But a single-word misspelling shares no
            # whole token with its target either ("Pomerarian"/"Pomeranian"),
            # so allow those through on a high character-level ratio alone.
            if not (len(a) == 1 and len(b) <= 2 and r >= 0.85):
                continue
        if a <= b or b <= a:
            r = max(r, 0.95)
        if r > score:
            best, score = c, r
    return (best, score) if score >= MATCH_MIN else (None, score)


def load_reported(path):
    out = {}
    for i, line in enumerate(open(path)):
        if i == 0:
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 4:
            out[f[0]] = (f[1], f[2], f[3])
    return out


def load_predictions(path, names):
    out = {}
    for i, line in enumerate(open(path)):
        if i == 0:
            continue
        f = line.rstrip('\n').split('\t')
        top = []
        for j in (1, 3, 5):
            if j < len(f) and f[j]:
                top.append((names.get(f[j], f[j]), float(f[j + 1] or 0)))
        out[f[0].split('-')[-1]] = top
    return out


def main():
    rep_path, pred_path, names_path = sys.argv[1:4]
    names = json.load(open(names_path))
    vocab = sorted(set(names.values()))
    reported = load_reported(rep_path)
    predicted = load_predictions(pred_path, names)

    pure, mixed, designer, unmapped = [], [], [], []
    for n in sorted(set(reported) & set(predicted), key=int):
        _, raw, sex = reported[n]
        top = predicted[n]
        key = norm(raw)

        cross = next((v for k, v in CROSSES.items() if k in key), None)
        if cross:
            designer.append((n, raw, list(cross), top))
            continue

        is_mix = bool(MIX_RE.search(raw))
        parts = [p for p in SPLIT_RE.split(MIX_RE.sub(' ', raw)) if norm(p)]
        exp = []
        for p in parts:
            c, _ = canonical(p, vocab)
            if c and c not in exp:
                exp.append(c)
        m = re.search(r'\(([^)]*)\)', raw)
        if m:
            c, _ = canonical(m.group(1), vocab)
            if c and c not in exp:
                exp.append(c)

        if not exp:
            unmapped.append((n, raw, top))
        elif is_mix or len(exp) > 1:
            mixed.append((n, raw, exp, top))
        else:
            pure.append((n, raw, exp, top))

    def names_of(top):
        return [p for p, _ in top]

    h1 = [r for r in pure if r[3] and r[2][0] == r[3][0][0]]
    h3 = [r for r in pure if r[2][0] in names_of(r[3])]

    print(f"dogs compared: {len(reported.keys() & predicted.keys())}\n")
    print(f"PUREBRED reports        n={len(pure)}")
    if pure:
        print(f"   top-1 correct      {len(h1):>3}/{len(pure)}  ({100*len(h1)/len(pure):.0f}%)")
        print(f"   reported in top-3  {len(h3):>3}/{len(pure)}  ({100*len(h3)/len(pure):.0f}%)")

    mall = [r for r in mixed if all(e in names_of(r[3]) for e in r[2])]
    many = [r for r in mixed if any(e in names_of(r[3]) for e in r[2])]
    print(f"\nMIXED / multi-breed reports  n={len(mixed)}")
    if mixed:
        print(f"   all named breeds in top-3  {len(mall):>3}/{len(mixed)}  ({100*len(mall)/len(mixed):.0f}%)")
        print(f"   >=1 named breed in top-3   {len(many):>3}/{len(mixed)}  ({100*len(many)/len(mixed):.0f}%)")

    dall = [r for r in designer if all(any(e.split()[0] in p for p in names_of(r[3])) for e in r[2])]
    dany = [r for r in designer if any(any(e.split()[0] in p for p in names_of(r[3])) for e in r[2])]
    print(f"\nDESIGNER CROSSES (no Parker label)  n={len(designer)}")
    if designer:
        print(f"   both parents in top-3      {len(dall):>3}/{len(designer)}")
        print(f"   >=1 parent in top-3        {len(dany):>3}/{len(designer)}")

    print(f"\n--- purebred misses ({len(pure)-len(h3)}) ---")
    for r in pure:
        if r not in h3:
            got = ', '.join(f'{p} {v:.2f}' for p, v in r[3][:2])
            print(f"  #{r[0]:<4} {r[1][:30]:<30} expected={r[2][0][:26]:<26} got={got}")

    print(f"\n--- unmatched to Parker vocabulary ({len(unmapped)}) ---")
    for n, raw, top in unmapped:
        print(f"  #{n:<4} {raw[:38]:<38} got={top[0][0]} {top[0][1]:.2f}" if top else f"  #{n} {raw}")


if __name__ == '__main__':
    main()
