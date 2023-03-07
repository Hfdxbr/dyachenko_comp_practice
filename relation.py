from pathlib import Path

import pandas as pd


def find_relation(pattern: str):
    dfs = {}
    for file in Path('.').glob(f'u_{pattern}*.csv'):
        N, M = map(lambda x: int(x[1:]), file.stem.split('_')[-2:])
        df = pd.read_csv(file)
        dfs[(N, M)] = df
        df.round(6).to_csv(file, index=False)
    keys = sorted(dfs.keys())
    out = []
    t_col = r'\(t\)'
    for (M1, N1), (M2, N2) in zip(keys, keys[1:]):
        print(M2, N2, '--', M1, N1)
        df = dfs[(M2, N2)] - dfs[(M1, N1)]
        out.append(df)
    assert len(out) == 2
    df = out[0] / out[1]
    df[t_col] = dfs[keys[0]][t_col]
    df.rename(columns=lambda x: x.replace(r'\left.u\right\vert_', ''), inplace=True)
    df.round(6).astype(str).to_csv(f'rel_{pattern}.csv', index=False)


find_relation('a0.1_b0.1')
find_relation('a1.0_b0.1')
find_relation('a1.0_b1.0')