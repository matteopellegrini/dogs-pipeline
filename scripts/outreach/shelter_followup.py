#!/usr/bin/env python3
"""Shelter outreach follow-up (one bump, same thread). Dry-run by default; --send sends.

Reads  ~/Downloads/shelter_outreach_list_2026-10-06.csv (+ sendlog 2026-10-06),
writes ~/Downloads/shelter_outreach_followup_sendlog_2026-10-19.json and fills followup_on.
"""
import csv, json, os, re, sys, time, urllib.request, urllib.error
from datetime import datetime, timezone

HOME = os.path.expanduser('~')
CSV = f'{HOME}/Downloads/shelter_outreach_list_2026-10-06.csv'
LOG0 = f'{HOME}/Downloads/shelter_outreach_sendlog_2026-10-06.json'
LOG1 = f'{HOME}/Downloads/shelter_outreach_followup_sendlog_2026-10-19.json'
API = 'https://my.prosperk9.com/api/admin/send-email'
KEY = sys.argv[sys.argv.index('--key') + 1] if '--key' in sys.argv else ''
SEND = '--send' in sys.argv
SUBJECT = "Re: A free place for your adopters to keep their dog's vet records"
TODAY = '2026-10-19'

BODY = """{greeting},

Bumping this once in case it got buried. I'd also like to add something concrete: I'll send you a free whole-genome DNA test for one dog of your choosing, ideally one who has been with you too long. You get the breed mix with chromosome-by-chromosome painting, 450 health conditions and drug sensitivities, and a one-page vet summary. Most rescues post the reveal, and it tends to get a dog noticed. A real report: https://my.prosperk9.com/demo

If you'd like the kit, reply with the dog's name and a shipping address and it goes out this week. The free Healthbook blurb for your adoption emails is yours either way. If it isn't a fit, no need to reply and I won't write again.

Matteo Pellegrini
Geneticist and founder, ProsperK9
prosperk9.com
"""

mask = lambda e: re.sub(r'([A-Za-z0-9._%+-]{2})[A-Za-z0-9._%+-]*@([A-Za-z0-9.-]+)', r'\1***@\2', e)

def api(method, path='', body=None):
    req = urllib.request.Request(API + path, data=json.dumps(body).encode() if body else None,
                                 headers={'x-api-key': KEY, 'Content-Type': 'application/json', 'User-Agent': 'Mozilla/5.0'}, method=method)
    try:
        with urllib.request.urlopen(req, timeout=30) as r: return r.status, json.loads(r.read())
    except urllib.error.HTTPError as e: return e.code, {'error': e.read()[:200].decode(errors='replace')}

rows = list(csv.DictReader(open(CSV)))
fields = list(rows[0].keys())
log0 = {e['to'].strip().lower(): e for e in json.load(open(LOG0))}
done = {e['to'] for e in json.load(open(LOG1))} if os.path.exists(LOG1) else set()

plan, skipped = [], []
for r in rows:
    email = r['email'].strip().lower()
    e0 = log0.get(email)
    why = None
    if r['status'] != 'sent': why = f'status={r["status"]}'
    elif (r.get('reply') or '').strip(): why = 'replied'
    elif (r.get('followup_on') or '').strip(): why = f'already followed up {r["followup_on"]}'
    elif not e0 or not e0.get('ok'): why = 'no ok send on 2026-10-06'
    elif email in done: why = 'in follow-up log'
    (skipped if why else plan).append((r, e0, why))

print(f'{len(plan)} to send, {len(skipped)} skipped')
for r, e0, why in skipped: print(f'  skip {r["org"]:35s} {mask(r["email"])}  ({why})')

out = json.load(open(LOG1)) if os.path.exists(LOG1) else []
for r, e0, _ in plan:
    greeting = (r.get('greeting') or 'Hi there').strip()
    text = BODY.format(greeting=greeting)
    # Thread under the original: look up its Message-ID from Resend (best effort).
    st, d = api('GET', f'?id={e0["id"]}') if KEY else (0, {})
    mid = d.get('message_id') if st == 200 else None
    print(f'  {"SEND" if SEND else "plan"} {r["org"]:35s} {mask(r["email"]):28s} {greeting:12s} thread={"yes" if mid else "subject-only"}')
    if not SEND: continue
    body = {'to': r['email'].strip(), 'subject': SUBJECT, 'text': text}
    if mid: body.update({'in_reply_to': mid, 'references': mid})
    st, d = api('POST', '', body)
    ok = st == 200 and d.get('ok')
    out.append({'org': r['org'], 'to': r['email'].strip(), 'sent_at': datetime.now(timezone.utc).isoformat(), 'ok': bool(ok),
                'id': d.get('id'), 'in_reply_to': mid, 'subject': SUBJECT, **({} if ok else {'error': d.get('error')})})
    json.dump(out, open(LOG1, 'w'), indent=1)
    if ok: r['followup_on'] = TODAY
    else: print('    FAILED', st, d)
    time.sleep(1.5)

if SEND:
    with open(CSV, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    print(f'sent {sum(1 for e in out if e["ok"])}/{len(out)}; log {LOG1}; CSV followup_on filled')
