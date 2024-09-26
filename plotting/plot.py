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
        title_text = f"{df.iloc[0, 0]} - {df.iloc[0, 1]}"
        fig = px.scatter(df, x=x_axis, y='p99', hover_data=hover_fields, color=color_by)

    else:
        df_pareto = pareto_front(df, x_axis)
        title_text = f"{df_pareto.iloc[0, 0]} - {df_pareto.iloc[0, 1]}"
        fig = px.scatter(df_pareto, x=x_axis, y='p99', hover_data=hover_fields, color=color_by)

    fig.update_layout(title=title_text, title_x=0.5)
    fig.update_layout(plot_bgcolor='white')
    fig.update_xaxes(
        mirror=False,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    )
    fig.update_yaxes(
        mirror=False,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    )
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