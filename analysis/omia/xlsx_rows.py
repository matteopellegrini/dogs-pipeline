#!/usr/bin/env python3
"""Minimal .xlsx reader (stdlib only): xlsx_rows.py FILE [SHEET_INDEX] -> TSV on stdout. Lists sheets if no index."""
import re, sys, zipfile
import xml.etree.ElementTree as ET
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
RID = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'


def col(ref):
    n = 0
    for ch in re.match(r'[A-Z]+', ref).group():
        n = n * 26 + ord(ch) - 64
    return n - 1


def main():
    z = zipfile.ZipFile(sys.argv[1])
    shared = []
    if 'xl/sharedStrings.xml' in z.namelist():
        for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', NS):
            shared.append(''.join(t.text or '' for t in si.iter('{%s}t' % NS['m'])))
    wb = ET.fromstring(z.read('xl/workbook.xml'))
    rels = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    target = {r.get('Id'): r.get('Target') for r in rels}
    sheets = [(s.get('name'), 'xl/' + target[s.get(RID)].lstrip('/').removeprefix('xl/')) for s in wb.find('m:sheets', NS)]
    if len(sys.argv) < 3:
        for i, (n, p) in enumerate(sheets):
            print(i, n, p)
        return
    root = ET.fromstring(z.read(sheets[int(sys.argv[2])][1]))
    for row in root.iter('{%s}row' % NS['m']):
        vals = {}
        for c in row.findall('m:c', NS):
            v = c.find('m:v', NS)
            t = c.get('t')
            if t == 's' and v is not None:
                x = shared[int(v.text)]
            elif t == 'inlineStr':
                x = ''.join(e.text or '' for e in c.iter('{%s}t' % NS['m']))
            else:
                x = v.text if v is not None else ''
            vals[col(c.get('r'))] = x.replace('\t', ' ').replace('\n', ' ')
        if vals:
            print('\t'.join(vals.get(i, '') for i in range(max(vals) + 1)))


if __name__ == '__main__':
    main()
