import plotly.express as px
import pandas as pd

def plot(df, 
         x_axis, 
         hover_fields = ['p99', 'hit_chance', 'period', 'utilization'], 
         color_by = None, 
         filter_fn = None, 
         allpoints = True):
    
    if filter_fn is not None:
        df = filter_fn(df)
    
    if allpoints:
        fig = px.scatter(df, x=x_axis, y='p99', hover_data=hover_fields, color=color_by)
        return fig
    else:
        df_pareto = pareto_front(df, x_axis)
        fig = px.scatter(df_pareto, x=x_axis, y='p99', hover_data=hover_fields, color=color_by)
    return fig
import pandas as pd

def pareto_front(df, col):
    df_sorted = df.sort_values(by=col, ascending=True)
    pareto_front = []
    min_dev = float('inf')
    for index, row in df_sorted.iterrows():
        if row['p99'] < min_dev:
            pareto_front.append(row)
            min_dev = row['p99']  
    return pd.DataFrame(pareto_front)

def readfile(DATA_PATH, JOB_ID, FILE_NUM):
    FILES = [f"{DATA_PATH}/{JOB_ID}/{i}.csv" for i in FILE_NUM]
    dataframes = []
    for i, file in zip(FILE_NUM, FILES):
        ref = pd.read_csv(file, names=None, index_col=None)
        dataframes.append(ref)
    df = pd.concat(dataframes)
    return df