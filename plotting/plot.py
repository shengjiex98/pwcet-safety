import plotly.express as px
import pandas as pd
from typing import Callable

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

def plot(
        df: pd.DataFrame,
        x: str,
        hover_fields: list[str] = ["p99", "hit_chance", "period", "utilization"],
        color_by: str = None,
        filter_fn: Callable[[pd.Series], bool] = None,
        allpoints: bool = True,
        highlight_min: bool = False
    ):
    if filter_fn is not None:
        df = df.loc[df.apply(filter_fn, axis=1)]
    if allpoints:
        title_text = f"{df.iloc[0, 0]} - {df.iloc[0, 1]}"
        fig = px.scatter(df, x=x, y="p99", hover_data=hover_fields, color=color_by, color_discrete_sequence=colors)
    else:
        df_pareto = df.loc[df.groupby(x)["p99"].idxmin()]
        min_row = df_pareto.loc[df_pareto['p99'].idxmin()]
        title_text = f"{df_pareto.iloc[0, 0]} - {df_pareto.iloc[0, 1]}"
        fig = px.line(
            df_pareto, x=x, y="p99", hover_data=hover_fields, color=color_by, markers=True, color_discrete_sequence=colors
        )
    fig.update_traces(marker=dict(size=8), line=dict(width=3))
    if highlight_min:
        fig.add_scatter(
            x=[min_row[x]],
            y=[min_row['p99']],
            mode='markers',
            marker=dict(size=12, color='red', symbol='circle'),
            name='P*'
            # name=f'Min Point (period={min_row['period']:.2f}, hit_chance={min_row['hit_chance']:.2f})'
        ).update_layout(
            legend=dict(
                x=0.99,
                y=0.99,
                xanchor='right',
                yanchor='top',
                bgcolor='rgba(255, 255, 255, 0.5)'  # semi-transparent background to improve visibility
            )
        )
    # fig.add_shape(
    #     type="rect",
    #     xref="paper", yref="paper",  # Reference the full plot area
    #     x0=0, y0=0, x1=1, y1=1,      # Cover the full plot area
    #     line=dict(color="black", width=1)  # Border color and width
    # )
    fig.update_layout(
        font=dict(size=20),
        template="plotly_white"
    ).update_xaxes(
        rangemode="tozero"
    ).update_yaxes(
        rangemode="tozero"
    )
    # fig.update_layout(title=title_text, title_x=0.5)
    
    # fig.update_layout(plot_bgcolor="white")
    # fig.update_xaxes(
    #     mirror=False, ticks="outside", showline=True, linecolor="black", gridcolor="lightgrey"
    # )
    # fig.update_yaxes(
    #     mirror=False, ticks="outside", showline=True, linecolor="black", gridcolor="lightgrey"
    # )
    return fig
