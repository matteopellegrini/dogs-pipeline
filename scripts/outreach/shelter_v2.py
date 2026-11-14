#!/usr/bin/env python3
"""Shelter outreach v2 (Healthbook only, LA rescues). Dry-run by default.

  shelter_v2.py --batch 1|2 [--send] [--key K]     first/second half of status=todo rows
  shelter_v2.py --followup [--send] [--key K]      day-6 bump to sent rows with empty reply and no followup_on
                                                    (only rows whose sent_on is >= 6 days ago)
List  ~/Downloads/shelter_outreach_v2_2026-10-19.csv (25 rows: 11 LA + 14 breed rescues)
Log   ~/Downloads/shelter_outreach_v2_sendlog_2026-10-19.json
"""
import csv, json, os, re, sys, time, urllib.request, urllib.error
from datetime import datetime, timezone, date, timedelta

HOME = os.path.expanduser('~')
CSV = f'{HOME}/Downloads/shelter_outreach_v2_2026-10-19.csv'
OLD = f'{HOME}/Downloads/shelter_outreach_list_2026-10-06.csv'
LOG = f'{HOME}/Downloads/shelter_outreach_v2_sendlog_2026-10-19.json'
API = 'https://my.prosperk9.com/api/admin/send-email'
KEY = sys.argv[sys.argv.index('--key') + 1] if '--key' in sys.argv else ''
SEND = '--send' in sys.argv
FOLLOWUP = '--followup' in sys.argv
BATCH = int(sys.argv[sys.argv.index('--batch') + 1]) if '--batch' in sys.argv else 0
MAX_PER_DAY = 15
SUBJECT = 'Free health record for every dog you adopt out'
TODAY = date.today().isoformat()

BODY = """{greeting},

I'm Matteo, a geneticist in Los Angeles and the founder of ProsperK9. One question: how often does an adopter write back asking you to resend the vaccine record or the spay certificate?

We built a free fix for that, the ProsperK9 Healthbook. Adopters photograph the adoption paperwork once, it's kept for life alongside future vet visits, vaccines, labs and photos, and an AI Health Advisor explains what's in the records and what to ask the vet (guidance to discuss with a vet, never a diagnosis). Free, no card, nothing to buy. Here's what adopters get: prosperk9.com/pages/healthbook

If {name} would like to offer it, I'll send you a one-line blurb and a small image with your name on it for the welcome email you already send. Nothing on it sells anything. To be transparent, ProsperK9 also sells a DNA test; some Healthbook users choose to buy it later, but adopters never need it, and the card doesn't mention it.

If it sounds useful, reply "yes" with your name as you'd like it shown and I'll send the files.

Thank you for the work you do,
Matteo Pellegrini
Geneticist and founder, ProsperK9
prosperk9.com
"""

FOLLOWUP_BODY = """{greeting}, bumping this once. The free Healthbook blurb for your adoption emails is ready whenever you want it; a "yes" is all I need. If it isn't a fit, no need to reply and I won't write again.

Matteo
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
old = {r['email'].strip().lower() for r in csv.DictReader(open(OLD))} | {r['org'].strip().lower() for r in csv.DictReader(open(OLD))}
log = json.load(open(LOG)) if os.path.exists(LOG) else []
sent_today = sum(1 for e in log if e.get('ok') and e['sent_at'][:10] == TODAY)
first_id = {e['to'].lower(): e['id'] for e in log if e.get('ok') and e.get('kind') == 'initial'}

plan = []
if FOLLOWUP:
    cutoff = (date.today() - timedelta(days=6)).isoformat()
    for r in rows:
        if r['status'] != 'sent' or (r.get('reply') or '').strip() or (r.get('followup_on') or '').strip(): continue
        if not r.get('sent_on') or r['sent_on'] > cutoff: print(f'  wait {r["org"]:38s} sent {r.get("sent_on")} (<6 days)'); continue
        if r['email'].strip().lower() not in first_id: print(f'  skip {r["org"]:38s} no ok initial send in log'); continue
        plan.append(r)
else:
    todo = [r for r in rows if r['status'] == 'todo' and r['email'].strip().lower() not in old and r['org'].strip().lower() not in old
            and not (r.get('sent_on') or '').strip() and r['email'].strip().lower() not in first_id]
    # Batch 1 = first half of what is still to do; batch 2 = everything left (both capped at MAX_PER_DAY).
    half = (len(todo) + 1) // 2
    plan = todo[:half] if BATCH == 1 else todo if BATCH == 2 else []
    if BATCH not in (1, 2): print('need --batch 1 or --batch 2'); sys.exit(1)
    plan = plan[:max(0, MAX_PER_DAY - sent_today)]

print(f'{"FOLLOWUP" if FOLLOWUP else f"batch {BATCH}"}: {len(plan)} to send (already sent today: {sent_today}); {"SENDING" if SEND else "dry run"}')
for r in plan:
    greeting = (r.get('greeting') or 'Hi there').strip()
    name = (r.get('ecard_display_name') or r['org']).strip()
    text = FOLLOWUP_BODY.format(greeting=greeting) if FOLLOWUP else BODY.format(greeting=greeting, name=name)
    subject = f'Re: {SUBJECT}' if FOLLOWUP else SUBJECT
    print(f'  {"SEND" if SEND else "plan"} {r["org"][:38]:38s} {mask(r["email"]):28s} {greeting:13s} {name}')
    if not SEND: continue
    body = {'to': r['email'].strip(), 'subject': subject, 'text': text}
    mid = None
    if FOLLOWUP and KEY:
        st, d = api('GET', f'?id={first_id[r["email"].strip().lower()]}')
        mid = d.get('message_id') if st == 200 else None
        if mid: body.update({'in_reply_to': mid, 'references': mid})
    st, d = api('POST', '', body)
    ok = st == 200 and bool(d.get('ok'))
    log.append({'org': r['org'], 'to': r['email'].strip(), 'sent_at': datetime.now(timezone.utc).isoformat(), 'ok': ok, 'id': d.get('id'),
                'subject': subject, 'kind': 'followup' if FOLLOWUP else 'initial', **({'in_reply_to': mid} if mid else {}), **({} if ok else {'error': d.get('error')})})
    json.dump(log, open(LOG, 'w'), indent=1)
    if ok:
        if FOLLOWUP: r['followup_on'] = TODAY
        else: r['status'] = 'sent'; r['sent_on'] = TODAY
    else: print('    FAILED', st, d)
    time.sleep(1.5)

if SEND and plan:
    with open(CSV, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    print(f'done: {sum(1 for e in log if e["ok"])} ok of {len(log)} in log; CSV updated')
