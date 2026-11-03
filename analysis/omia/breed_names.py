"""Map free-text breed names (OMIA, commercial panels) onto breed panel codes.

match(name) returns the set of panel codes a published breed name refers to: an exact match on the normalised
name, otherwise every panel breed whose name contains it as whole words ("Poodle" -> Standard/Miniature/Toy,
"Dachshund" -> every Dachshund), otherwise the panel breed whose name it contains. ALIASES covers the names that
differ in wording. "Mixed breed" and similar return MIXED (the caller treats it as "seen outside purebreds").
"""
import re

MIXED = 'MIXED'
ALIASES = {
    'german shepherd dog': 'german shepherd', 'doberman pinscher': 'doberman', 'dobermann': 'doberman',
    'alsatian': 'german shepherd', 'bull mastiff': 'bullmastiff', 'shar pei': 'chinese shar pei',
    'chinese sharpei': 'chinese shar pei', 'sharpei': 'chinese shar pei', 'akita inu': 'akita',
    'japanese akita': 'akita', 'old english sheepdog': 'old english sheepdog', 'yorkie': 'yorkshire terrier',
    'westie': 'west highland white terrier', 'staffie': 'staffordshire bull terrier',
    'english bulldog': 'bulldog', 'british bulldog': 'bulldog', 'pit bull': 'american pit bull terrier',
    'pitbull': 'american pit bull terrier', 'cavalier': 'cavalier king charles spaniel',
    'king charles spaniel': 'cavalier king charles spaniel', 'cocker spaniel': 'cocker spaniel',
    'american cocker spaniel': 'cocker spaniel', 'springer spaniel': 'english springer spaniel',
    'wirehaired dachshund': 'dachshund', 'longhaired dachshund': 'dachshund', 'smooth dachshund': 'dachshund',
    'standard schnauzer': 'standard schnauzer', 'schnauzer': 'schnauzer', 'mexican hairless': 'xoloitzcuintli',
    'xolo': 'xoloitzcuintli', 'rough collie': 'collie', 'smooth collie': 'collie', 'sheltie': 'shetland sheepdog',
    'belgian shepherd': 'belgian', 'belgian shepherd dog': 'belgian', 'malinois': 'belgian malinois',
    'tervuren': 'belgian tervuren', 'groenendael': 'belgian sheepdog', 'labrador': 'labrador retriever',
    'golden': 'golden retriever', 'toller': 'nova scotia duck tolling retriever', 'dogue de bordeaux': 'dogue de bordeaux',
    'french mastiff': 'dogue de bordeaux', 'great pyrenees': 'great pyrenees', 'pyrenean mountain dog': 'great pyrenees',
    'bernese': 'bernese mountain dog', 'bernese mountain dogs': 'bernese mountain dog', 'parson russell': 'parson russell terrier',
    'jack russell': 'jack russell terrier', 'american staffordshire': 'american staffordshire terrier',
    'wire fox terrier': 'fox terrier', 'smooth fox terrier': 'fox terrier', 'miniature poodle': 'miniature poodle',
    'st bernard': 'saint bernard', 'irish red and white setter': 'irish setter',
    'welsh corgi cardigan': 'cardigan welsh corgi', 'welsh corgi pembroke': 'pembroke welsh corgi',
    'cairn terrer': 'cairn terrier', 'italian spinone': 'spinone italiano',
    # several panel breeds (proxies where the panel lacks the breed itself)
    'german pointer': ['german shorthaired pointer', 'german wirehaired pointer'],
    'russell type terriers': ['jack russell terrier', 'parson russell terrier'],
    'labradoodle': ['labrador retriever', 'standard poodle', 'miniature poodle', 'toy poodle'],
    'lapponian herder': ['finnish lapphund', 'swedish lapphund'],
    'alaskan husky': ['siberian husky'],
}
MIXED_WORDS = ('mixed breed', 'mixed', 'crossbred', 'mongrel', 'multiple breeds', 'various', 'many breeds', 'many',
               'multi breed', 'all dogs', 'western breeds', 'eastern breeds')
IGNORE = ('et al', 'etc', 'dogs')  # list fragments, not breeds

def norm(s):
    s = re.sub(r'\(dog\)|\(canis[^)]*\)', ' ', s.lower())
    s = re.sub(r"[^a-z ]", ' ', s.replace('-', ' '))
    return re.sub(r'\s+', ' ', s).strip()

class BreedMatcher:
    def __init__(self, names):
        """names: {code: display name}."""
        self.codes = {c: norm(n) for c, n in names.items()}
        self.by_name = {}
        for c, n in self.codes.items():
            self.by_name.setdefault(n, set()).add(c)

    def match(self, name):
        q = norm(name)
        if not q:
            return set()
        if q in IGNORE or q.split(' ')[-1] == 'dogs' and q[0].isdigit():
            return set()
        if any(q == w or q.startswith(w + ' ') for w in MIXED_WORDS):
            return {MIXED}
        q = ALIASES.get(q, q)
        if isinstance(q, list):
            return set().union(*(self.match(x) for x in q))
        if q in self.by_name:
            return set(self.by_name[q])
        pat = re.compile(r'(^| )' + re.escape(q) + r'( |$)')
        hits = {c for c, n in self.codes.items() if pat.search(n)}
        if hits:
            return hits
        hits = {c for c, n in self.codes.items() if len(n) > 3 and re.search(r'(^| )' + re.escape(n) + r'( |$)', q)}
        if hits:
            longest = max(len(self.codes[c]) for c in hits)
            return {c for c in hits if len(self.codes[c]) == longest}
        if q.endswith('s') and q[:-1] in self.by_name:
            return set(self.by_name[q[:-1]])
        return set()
