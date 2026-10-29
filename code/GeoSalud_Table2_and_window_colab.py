# -*- coding: utf-8 -*-
"""
GeoSalud Cotopaxi - Table 2 (raw series) and eruptive-period window (weeks 40-52).
Run in Google Colab or locally.  INPUT: weekly_panel.csv  (columns: canton, year, week, so2_raw, discharges)
OUTPUT: Table2_reproduced.csv, Eruptive_window_40_52.csv  (printed as well)
Definitions (as in the Methods): raw, unreplaced SO2; weekly pairs with a missing value are dropped;
Pearson and Spearman correlations; Benjamini-Hochberg FDR (q = 0.05) across the seven cantons within each year;
Fisher r-to-z test for independent correlations (2022 vs 2024); eruptive window = ISO weeks 40-52 (13 weeks).
"""
import os, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

URL = 'https://raw.githubusercontent.com/Santiago031289/cotopaxi-so2/main/data/weekly_panel.csv'
df = pd.read_csv('weekly_panel.csv') if os.path.exists('weekly_panel.csv') else pd.read_csv(URL)
ORDER = ['Latacunga', 'Pujilí', 'Salcedo', 'La Maná', 'Saquisilí', 'Pangua', 'Sigchos']
CORE = ['Latacunga', 'Pujilí', 'Saquisilí', 'Salcedo']

def corr(g):
    h = g.dropna(subset=['so2_raw', 'discharges'])
    r, p = stats.pearsonr(h.so2_raw, h.discharges); rs, ps = stats.spearmanr(h.so2_raw, h.discharges)
    return len(h), r, p, rs, ps

rows = []
for c in ORDER:
    for y in (2022, 2024):
        n, r, p, rs, ps = corr(df[(df.canton == c) & (df.year == y)])
        rows.append(dict(canton=c, year=y, n=n, pearson_r=r, pearson_p=p, spearman_r=rs, spearman_p=ps))
T2 = pd.DataFrame(rows)
for y in (2022, 2024):                                   # BH-FDR within each year (seven cantons)
    m = T2.year == y
    T2.loc[m, 'pearson_p_BH'] = multipletests(T2.loc[m, 'pearson_p'], method='fdr_bh')[1]
    T2.loc[m, 'significant_after_FDR'] = T2.loc[m, 'pearson_p_BH'] < 0.05
def fisher_p(r1, n1, r2, n2):
    z = (np.arctanh(r1) - np.arctanh(r2)) / np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3)); return 2 * (1 - stats.norm.cdf(abs(z)))
fp = {}
for c in ORDER:
    a = T2[(T2.canton == c) & (T2.year == 2022)].iloc[0]; b = T2[(T2.canton == c) & (T2.year == 2024)].iloc[0]
    fp[c] = fisher_p(a.pearson_r, a.n, b.pearson_r, b.n)
T2['fisher_p_2022_vs_2024'] = T2.canton.map(fp)
T2 = T2.round(4); T2.to_csv('Table2_reproduced.csv', index=False)

rows = []
for c in CORE:
    for y in (2022, 2024):
        g = df[(df.canton == c) & (df.year == y) & df.week.between(40, 52)]
        n, r, p, rs, ps = corr(g)
        rows.append(dict(canton=c, year=y, n_weeks=n, pearson_r=r, pearson_p=p, discharges_w40_52=int(g.discharges.sum())))
W = pd.DataFrame(rows).round(4); W.to_csv('Eruptive_window_40_52.csv', index=False)

pd.set_option('display.width', 220)
print('=== Table 2 (raw series) ==='); print(T2.to_string(index=False))
print('\n=== Eruptive window, weeks 40-52 ==='); print(W.to_string(index=False))
