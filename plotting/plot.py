import plotly.express as px
import pandas as pd
import os

def plot(df, 
         x_axis, 
         hover_fields = ['dev', 'q', 'period', 'utilization'], 
         color_by = None, 
         filter_fn = None, 
         allpoints = True):
    if filter_fn is not None:
        df = df[filter_fn(df)]
    

    if allpoints:
        fig = px.scatter(df, x=x_axis, y='dev', hover_data=hover_fields, color=color_by)
        return fig
    else:
        df_sorted = df.sort_values(by=df[x_axis])
        pareto_values = pareto_front(df_sorted, df_sorted['dev'])
        df_pareto = df_sorted[df_sorted['dev'].isin(pareto_values)]
        fig = px.scatter(df_pareto, x=x_axis, y='dev', hover_data=hover_fields, color=color_by)
    return fig

def pareto_front(df, col):
    pareto_front = []
    min = float('-inf')  
    for value in df[col]:
        if value <= min:
            pareto_front.append(value)
            min = value
    return pareto_front