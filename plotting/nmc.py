import plotly.express as px
import pandas as pd
import os

JOB_ID = 46080931
DATA_PATH = os.path.join(os.path.dirname(__file__), f"../data-proxy/nmc-grid/{JOB_ID}/")
FILE_NUM = range(1, 100+1)

FILES = [f"{DATA_PATH}/{JOB_ID}-{i}.csv" for i in FILE_NUM]

COL_NAMES = ['quantile', 'period', 'p99_low', 'p99', 'p99_high']

dataframes = []
for i, file in zip(FILE_NUM, FILES):
    ref = pd.read_csv(file, header=None, names=COL_NAMES, index_col=None)
    ref['util'] = range(100-len(ref)+1, 101)
    dataframes.append(ref)
df = pd.concat(dataframes)

fig_util = px.line(df.groupby('util')['p99'].min().reset_index(), x='util', y='p99', markers=True)
fig_util.show()

fig_ref = px.line(df.loc[(df['p99'] < 2) & (df['util'] == 100)], x='period', y='p99', hover_data=['quantile', 'period', 'p99'], markers=True)
fig_ref.show()
