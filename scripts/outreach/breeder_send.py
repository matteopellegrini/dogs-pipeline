#!/usr/bin/env python3
"""Breeder outreach (Golden Retriever + French Bulldog), google_ads BREEDER_SEND_REQUEST.md. Dry-run by default.

  breeder_send.py --day 1|2|3 [--send] [--key K]   initial email: day 1 = 10 rows, days 2-3 = 9, both breeds mixed
  breeder_send.py --followup [--send] [--key K]    day-6 bump, same subject thread, rows with empty reply only
List ~/Downloads/breeder_outreach_2026-10-27.csv (status=ok rows); log ~/Downloads/breeder_outreach_sendlog_2026-10-27.json.
Claims follow google_ads/BREEDER_PART_A_ANSWERS.md: only panel-imputed variants are named; COI yes; no parentage;
no OFA/AKC/club acceptance; no "only whole-genome".
"""
import csv, json, os, re, sys, time, urllib.request, urllib.error, urllib.parse
from datetime import datetime, timezone, date, timedelta

HOME = os.path.expanduser('~')
CSV = f'{HOME}/Downloads/breeder_outreach_2026-10-27.csv'
LOG = f'{HOME}/Downloads/breeder_outreach_sendlog_2026-10-27.json'
API = 'https://my.prosperk9.com/api/admin'
args = sys.argv[1:]
KEY = args[args.index('--key') + 1] if '--key' in args else ''
SEND = '--send' in args; FOLLOWUP = '--followup' in args
DAY = int(args[args.index('--day') + 1]) if '--day' in args else 0
PER_DAY = {1: 10, 2: 9, 3: 9}
TODAY = date.today().isoformat()

VARIANTS = {'Golden Retriever': 'prcd-PRA and degenerative myelopathy',
            'French Bulldog': 'hyperuricosuria and degenerative myelopathy'}
SUBJECT = "Whole-genome testing for {kennel_pos} next litter"
possessive = lambda k: k + ("'" if k.endswith("s") else "'s")
BODY = """{greeting},

I'm Matteo Pellegrini, a geneticist and the founder of ProsperK9. I saw that {kennel} health-tests its dogs, so I wanted to tell you about what we do differently: we sequence a dog's whole genome from one cheek swab instead of reading a marker chip.

For a {breed} that means about 450 inherited conditions and drug sensitivities, including {variants}, a genome-wide inbreeding coefficient, an oral microbiome report, and a one-page summary for the vet. The genome is re-analysed free whenever the science moves, so a test you run today keeps getting more useful.

For breeders the price is $125 a kit, or five for $625, shipping included. Your puppy buyers can also start a free ProsperK9 Healthbook, a lifelong record for vet visits, vaccines and notes.

A real report, so you can judge it yourself: https://my.prosperk9.com/demo?src=breeder_email

If it would help, I'm happy to send one kit for one of your own dogs. Just reply with the dog's name.

Matteo Pellegrini
Geneticist and founder, ProsperK9
prosperk9.com
ProsperK9, Building 114 CNSI, Los Angeles, CA 90095
"""
FOLLOWUP_BODY = """{greeting}, bumping this once. The offer of a kit for one of your own dogs stands if you'd like to see a whole-genome report first. If it isn't a fit, no need to reply and I won't write again.

Matteo
"""
mask = lambda e: re.sub(r'([A-Za-z0-9._%+-]{2})[A-Za-z0-9._%+-]*@([A-Za-z0-9.-]+)', r'\1***@\2', e)

def greeting(name):
    n = (name or '').strip()
    if not n or re.search(r'&|\band\b|,|/', n): return 'Hi there'
    first = [w for w in n.split() if not re.match(r'^(dr|mr|mrs|ms)\.?$', w, re.I)]
    return f'Hi {first[0]}' if first and len(first[0]) > 2 and not first[0].endswith('.') else 'Hi there'

def api(method, path, body=None):
    req = urllib.request.Request(API + path, data=json.dumps(body).encode() if body else None, method=method,
                                 headers={'x-api-key': KEY, 'Content-Type': 'application/json', 'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req, timeout=30) as r: return r.status, json.loads(r.read())
    except urllib.error.HTTPError as e: return e.code, {'error': e.read()[:200].decode(errors='replace')}

rows = list(csv.DictReader(open(CSV))); fields = list(rows[0].keys())
log = json.load(open(LOG)) if os.path.exists(LOG) else []
sent_to = {e['to'].lower() for e in log if e.get('ok') and e.get('kind') == 'initial'}
followed = {e['to'].lower() for e in log if e.get('ok') and e.get('kind') == 'followup'}

if FOLLOWUP:
    cutoff = (date.today() - timedelta(days=6)).isoformat()
    plan = [r for r in rows if r['status'] == 'ok' and r['sent'] and r['sent'] <= cutoff and not r['reply'].strip()
            and r['email'].lower() not in followed]
else:
    if DAY not in PER_DAY: print('need --day 1|2|3 or --followup'); sys.exit(1)
    todo = [r for r in rows if r['status'] == 'ok' and not r['sent'] and r['email'].lower() not in sent_to]
    g = [r for r in todo if r['breed'] == 'Golden Retriever']; f = [r for r in todo if r['breed'] == 'French Bulldog']
    n = PER_DAY[DAY]; nf = min(len(f), round(n * len(f) / max(1, len(todo)))) if DAY < 3 else len(f)
    plan = (f[:nf] + g[:n - nf])[:n] if DAY < 3 else (f + g)[:15]

print(f'{"FOLLOWUP" if FOLLOWUP else f"day {DAY}"}: {len(plan)} rows ({sum(r["breed"]=="Golden Retriever" for r in plan)} Golden, {sum(r["breed"]=="French Bulldog" for r in plan)} French Bulldog); {"SENDING" if SEND else "dry run"}')
for r in plan:
    gr = greeting(r['breeder_name']); subj = SUBJECT.format(kennel_pos=possessive(r['kennel']))
    if FOLLOWUP: subj, text = 'Re: ' + subj, FOLLOWUP_BODY.format(greeting=gr)
    else: text = BODY.format(greeting=gr, kennel=r['kennel'], breed=r['breed'], variants=VARIANTS[r['breed']])
    print(f'  {"SEND" if SEND else "plan"} {r["breed"][:6]:6s} {r["kennel"][:34]:34s} {mask(r["email"]):30s} {gr}')
    if not SEND: continue
    st, d = api('GET', '/optout?' + urllib.parse.urlencode({'email': r['email']}))
    if st == 200 and (d.get('opted_out') or d.get('optedOut')):
        print('    SKIP opted out'); continue
    st, d = api('POST', '/send-email', {'to': r['email'].strip(), 'subject': subj, 'text': text})
    ok = st == 200 and bool(d.get('ok'))
    log.append({'kennel': r['kennel'], 'to': r['email'].strip(), 'sent_at': datetime.now(timezone.utc).isoformat(), 'ok': ok,
                'id': d.get('id'), 'subject': subj, 'kind': 'followup' if FOLLOWUP else 'initial', **({} if ok else {'error': d.get('error')})})
    json.dump(log, open(LOG, 'w'), indent=1)
    if ok:
        if FOLLOWUP: r['followup_on'] = TODAY
        else: r['sent'] = TODAY
    else: print('    FAILED', st, d)
    time.sleep(1.5)
if SEND and plan:
    with open(CSV, 'w', newline='') as fh: w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)
    print(f'done: {sum(1 for e in log if e["ok"])} ok of {len(log)} in log; CSV updated')
if not plan: pass
elif not SEND:
    r = plan[0]; print('\n--- sample email ---\nSubject:', SUBJECT.format(kennel_pos=possessive(r['kennel'])) if not FOLLOWUP else 'Re: ...')
    print(BODY.format(greeting=greeting(r['breeder_name']), kennel=r['kennel'], breed=r['breed'], variants=VARIANTS[r['breed']]) if not FOLLOWUP else FOLLOWUP_BODY.format(greeting=greeting(r['breeder_name'])))
