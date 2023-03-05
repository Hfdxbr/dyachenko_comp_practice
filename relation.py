from pathlib import Path

import pandas as pd


def find_relation(pattern: str):
    dfs = {}
    for file in Path('.').glob(f'u_{pattern}*.csv'):
        N, M = map(lambda x: int(x[1:]), file.stem.split('_')[-2:])
        dfs[(N, M)] = pd.read_csv(file)
    keys = sorted(dfs.keys())
    out = []
    t_col = r'\(t\)'
    for (M1, N1), (M2, N2) in zip(keys, keys[1:]):
        print(M2, N2, '--', M1, N1)
        df = (dfs[(M2, N2)] / dfs[(M1, N1)]).round(4)
        df[t_col] = dfs[(M1, N1)][t_col]
        out.append(df)
    for df in out[1:]:
        for col in df.columns:
            if col == t_col: continue
            out[0][col] = out[0][col].astype(str) + ';' + df[col].astype(str)
    out[0].rename(columns=lambda x: x.replace('u', r"\frac{u_{MN}}{u_{M'N'}}"))
    out[0].to_csv(f'rel_{pattern}.csv', index=False)


find_relation('a0.1_b1.0')
find_relation('a1.0_b0.1')