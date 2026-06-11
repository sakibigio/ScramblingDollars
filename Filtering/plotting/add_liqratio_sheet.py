#!/usr/bin/env python3
"""
Add a 'LiquidityRatios' worksheet to data/LFX_datainputs.xlsx WITHOUT disturbing
any existing part (worksheets, 12 chartsheets, 66 charts, drawings, styles).

Done by surgical zip-level insertion: every original entry is copied byte-for-byte;
only [Content_Types].xml, xl/workbook.xml and xl/_rels/workbook.xml.rels are amended,
and one new part xl/worksheets/sheet24.xml is appended.

Series written (all in percent):
  FEDS-Note HQLA/total-assets, Standard/Modified/Non-LCR banks  (linear & pchip extended)
  Model e^mu_us, e^mu_eu  (empirical liquidity ratio from LFX_data.mat)
Row 1 = provenance/construction reference.
"""
import os, re, zipfile, shutil
import numpy as np
import scipy.io
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
XLSX = os.path.join(DATA, "LFX_datainputs.xlsx")

SHEET_NAME = "LiquidityRatios"
NEW_SHEET_FILE = "sheet24.xml"
NEW_RID = "rId46"
NEW_SHEETID = 41

PROVENANCE = (
    "Liquidity-ratio series (all in PERCENT, HQLA/total-assets). "
    "FEDS-Note columns: hand-digitized from FRB FEDS Note (Feb 2020) "
    "'The Liquidity Coverage Ratio and Corporate Liquidity Management', "
    "Standard/Modified/Non-LCR bank groups, quarterly 2003Q1-2018Q3 (values read "
    "off the figure, accuracy ~+/-0.5pp), interpolated to monthly by linear and pchip, "
    "padded flat to the 2001-01..2025-12 grid (padded regions carry no information "
    "beyond the endpoint level). Model e^mu_us / e^mu_eu: empirical liquidity ratio "
    "ingested by the filter, exp(mu) from data/LFX_data.mat (see load_data.m), "
    "monthly 2001-01..2025-07. Built 2026-06-04 via plotting/add_liqratio_sheet.py; "
    "sources hqla_assets_fedsnote_monthly_{linear,pchip}_extended.csv & LFX_data.mat."
)

def read_feds(fname):
    out = {}
    with open(os.path.join(DATA, fname)) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Date"):
                continue
            d, a, b, c = line.split(",")
            out[d] = (float(a), float(b), float(c))
    return out

lin = read_feds("hqla_assets_fedsnote_monthly_linear_extended.csv")
pch = read_feds("hqla_assets_fedsnote_monthly_pchip_extended.csv")

mat = scipy.io.loadmat(os.path.join(DATA, "LFX_data.mat"))
mu_us = np.asarray(mat["mu_us"]).ravel()
mu_eu = np.asarray(mat["mu_eu"]).ravel()
mu = {}
for k in range(len(mu_us)):
    d = "%04d-%02d" % (2001 + k // 12, k % 12 + 1)
    mu[d] = (float(np.exp(mu_us[k]) * 100.0), float(np.exp(mu_eu[k]) * 100.0))

# union grid = the extended FEDS grid (2001-01..2025-12), sorted
grid = sorted(lin.keys())

HEADERS = ["Date",
           "FEDS_Standard_LCR_linear", "FEDS_Modified_LCR_linear", "FEDS_Non_LCR_linear",
           "FEDS_Standard_LCR_pchip", "FEDS_Modified_LCR_pchip", "FEDS_Non_LCR_pchip",
           "Model_exp_mu_US", "Model_exp_mu_EU"]

def col(i):  # 0-based -> A,B,...
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s

def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))

def inline(cellref, text):
    return '<c r="%s" t="inlineStr"><is><t xml:space="preserve">%s</t></is></c>' % (cellref, esc(text))

def num(cellref, val):
    return '<c r="%s"><v>%s</v></c>' % (cellref, repr(round(val, 6)))

rows = []
# row 1: provenance (single cell A1)
rows.append('<row r="1">%s</row>' % inline("A1", PROVENANCE))
# row 2: headers
hc = "".join(inline(col(j) + "2", h) for j, h in enumerate(HEADERS))
rows.append('<row r="2">%s</row>' % hc)
# data rows from row 3
for ri, d in enumerate(grid):
    r = ri + 3
    cells = [inline("A%d" % r, d)]
    a, b, c = lin[d]
    cells += [num("B%d" % r, a), num("C%d" % r, b), num("D%d" % r, c)]
    a, b, c = pch[d]
    cells += [num("E%d" % r, a), num("F%d" % r, b), num("G%d" % r, c)]
    if d in mu:
        us, eu = mu[d]
        cells += [num("H%d" % r, us), num("I%d" % r, eu)]
    rows.append('<row r="%d">%s</row>' % (r, "".join(cells)))

sheet_xml = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
    'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
    '<sheetData>' + "".join(rows) + '</sheetData></worksheet>'
)

# --- repackage zip, amending only the three index parts ---------------------
TMP = XLSX + ".tmp"
with zipfile.ZipFile(XLSX, "r") as zin, \
     zipfile.ZipFile(TMP, "w", zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == "[Content_Types].xml":
            ins = ('<Override PartName="/xl/worksheets/%s" '
                   'ContentType="application/vnd.openxmlformats-officedocument.'
                   'spreadsheetml.worksheet+xml"/>' % NEW_SHEET_FILE).encode()
            data = data.replace(b"</Types>", ins + b"</Types>")
        elif item.filename == "xl/_rels/workbook.xml.rels":
            ins = ('<Relationship Id="%s" Type="http://schemas.openxmlformats.org/'
                   'officeDocument/2006/relationships/worksheet" '
                   'Target="worksheets/%s"/>' % (NEW_RID, NEW_SHEET_FILE)).encode()
            data = data.replace(b"</Relationships>", ins + b"</Relationships>")
        elif item.filename == "xl/workbook.xml":
            ins = ('<sheet name="%s" sheetId="%d" r:id="%s"/>'
                   % (SHEET_NAME, NEW_SHEETID, NEW_RID)).encode()
            data = data.replace(b"</sheets>", ins + b"</sheets>")
        zout.writestr(item, data)
    zout.writestr("xl/worksheets/" + NEW_SHEET_FILE, sheet_xml)

shutil.move(TMP, XLSX)
print("Added sheet '%s' (%s) with %d data rows (grid %s..%s)."
      % (SHEET_NAME, NEW_SHEET_FILE, len(grid), grid[0], grid[-1]))
