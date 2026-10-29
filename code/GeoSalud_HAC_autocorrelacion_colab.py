# GeoSalud - Supplementary Table S8: GLMs with autocorrelation-robust (Newey-West HAC) standard errors
# Pegar en una celda de Colab. Lee el panel semanal real del repositorio (o weekly_panel.csv local si ya lo subiste a Colab).
import os, numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats

URL = 'https://raw.githubusercontent.com/Santiago031289/cotopaxi-so2/main/data/weekly_panel.csv'
df = pd.read_csv('weekly_panel.csv') if os.path.exists('weekly_panel.csv') else pd.read_csv(URL)
CORE = ['Latacunga', 'Pujilí', 'Saquisilí', 'Salcedo']
MAXLAGS = 4   # ancho de ventana Newey-West (4 semanas); cambiar a 3 o 6 para ver sensibilidad

def run(series_col='so2_raw'):
    out = []
    for c in CORE:
        for y in (2022, 2024):
            g = df[(df.canton == c) & (df.year == y)].sort_values('week').reset_index(drop=True)
            s = g[series_col]
            lags = {f'lag{k}': s.shift(k) for k in range(4)}
            lags['lag_avg03'] = s.rolling(4, min_periods=4).mean()
            for name, v in lags.items():
                d = pd.DataFrame({'y': g.discharges, 'x': v, 'w': g.week.astype(float)}).dropna()
                d['w2'] = d.w ** 2
                iqr = d.x.quantile(.75) - d.x.quantile(.25)
                X = sm.add_constant(d[['x', 'w', 'w2']])
                fam = sm.families.Poisson()
                m = sm.GLM(d.y, X, family=fam).fit()
                disp = max(m.pearson_chi2 / m.df_resid, 1.0)
                b = m.params['x']
                se_qp = m.bse['x'] * np.sqrt(disp)
                h = sm.GLM(d.y, X, family=fam).fit(cov_type='HAC', cov_kwds={'maxlags': MAXLAGS})
                se_h = h.bse['x']
                # sensibilidad 2: se agrega el desenlace de la semana previa (término AR1) con SE quasi-Poisson
                d2 = d.copy(); d2['ylag'] = np.log1p(g.discharges.shift(1).reindex(d.index)); d2 = d2.dropna()
                X2 = sm.add_constant(d2[['x', 'w', 'w2', 'ylag']])
                m2 = sm.GLM(d2.y, X2, family=fam).fit(); disp2 = max(m2.pearson_chi2 / m2.df_resid, 1.0)
                b2 = m2.params['x']; se2 = m2.bse['x'] * np.sqrt(disp2)
                P = lambda bb, ss: 2 * (1 - stats.norm.cdf(abs(bb / ss)))
                out.append(dict(canton=c, year=y, lag=name, n=len(d), IQR=round(iqr, 6),
                                ER_pct=round((np.exp(b * iqr) - 1) * 100, 1),
                                HAC_CI_low=round((np.exp((b - 1.96 * se_h) * iqr) - 1) * 100, 1),
                                HAC_CI_high=round((np.exp((b + 1.96 * se_h) * iqr) - 1) * 100, 1),
                                p_quasi=round(P(b, se_qp), 4), p_HAC=round(P(b, se_h), 4),
                                ER_pct_AR1=round((np.exp(b2 * iqr) - 1) * 100, 1), p_AR1=round(P(b2, se2), 4)))
    return pd.DataFrame(out)

S8 = run('so2_raw')
S8.to_csv('Table_S8_HAC.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
print(S8[S8.year == 2022].to_string(index=False))
print('\n--- 2024: solo estimaciones con algun p<0.05 ---')
print(S8[(S8.year == 2024) & ((S8.p_quasi < .05) | (S8.p_HAC < .05) | (S8.p_AR1 < .05))].to_string(index=False))
print('\nSignificativos (p<0.05) en 2022 -> quasi-Poisson:', int((S8[S8.year == 2022].p_quasi < .05).sum()),
      '| HAC:', int((S8[S8.year == 2022].p_HAC < .05).sum()), '| AR1:', int((S8[S8.year == 2022].p_AR1 < .05).sum()))
