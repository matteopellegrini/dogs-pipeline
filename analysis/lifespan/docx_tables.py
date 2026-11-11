"""Dump the tables in a .docx as TSV (stdlib only): python3 -I docx_tables.py file.docx"""
import re, sys, zipfile
import xml.etree.ElementTree as ET
W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
root = ET.fromstring(zipfile.ZipFile(sys.argv[1]).read('word/document.xml'))
for ti, tbl in enumerate(root.iter(W + 'tbl')):
    print(f'## table {ti}')
    for tr in tbl.iter(W + 'tr'):
        cells = [''.join(t.text or '' for t in tc.iter(W + 't')).strip() for tc in tr.iter(W + 'tc')]
        print('\t'.join(re.sub(r'\s+', ' ', c) for c in cells))
