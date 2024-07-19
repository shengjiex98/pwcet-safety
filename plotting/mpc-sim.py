import plotly.express as px
import pandas as pd
import os

DATA_PATH = os.path.join(os.path.dirname(__file__), "../data-proxy/mpc-flags-uniform-period/")
FILE_NUM = range(1, 100+1)

REF_FILES = [f"{DATA_PATH}ref/43775543-{i}.csv" for i in FILE_NUM]
Y_FILES = [f"{DATA_PATH}y/43253817-{i}.csv" for i in FILE_NUM]

COL_NAMES = ['quantile', 'period', 'p99_low', 'p99', 'p99_high']

dfs = []
for i, file in zip(FILE_NUM, REF_FILES):
    ref = pd.read_csv(file, header=None, names=COL_NAMES, index_col=None)
    ref.insert(0, 'file_id', i)
    dfs.append(ref)
ref_df = pd.concat(dfs)

dfs = []
for i, file in zip(FILE_NUM, Y_FILES):
    y = pd.read_csv(file, header=None, names=COL_NAMES)
    y.insert(0, 'file_id', i)
    dfs.append(y)
y_df = pd.concat(dfs)

fig_ref = px.line(ref_df.loc[(ref_df['file_id'].isin([1, 8, 30])) & (ref_df['p99'] < 2500)], x='period', y='p99', color='file_id', hover_data=['quantile', 'period', 'p99'], title="Reference", markers=True)
fig_ref.show()

fig_y = px.line(y_df.loc[(y_df['file_id'].isin([1, 8, 30])) & (y_df['p99'] < 2500)], x='period', y='p99', color='file_id', hover_data=['quantile', 'period', 'p99'], title="Y", markers=True)
fig_y.show()
