import plotly.express as px
import pandas as pd
from typing import Callable


def plot(
    df: pd.DataFrame,
    x_axis: str,
    hover_fields: list[str] = ["p99", "hit_chance", "period", "utilization"],
    color_by: str = None,
    filter_fn: Callable[[pd.Series], bool] = None,
    allpoints: bool = True,
):
    if filter_fn is not None:
        df = df.loc[df.apply(filter_fn, axis=1)]
    if allpoints:
        title_text = f"{df.iloc[0, 0]} - {df.iloc[0, 1]}"
        fig = px.scatter(df, x=x_axis, y="p99", hover_data=hover_fields, color=color_by)
    else:
        df_pareto = df.loc[df.groupby(x_axis)["p99"].idxmin()]
        # df_pareto = pareto_front(df, x_axis)
        title_text = f"{df_pareto.iloc[0, 0]} - {df_pareto.iloc[0, 1]}"
        fig = px.line(
            df_pareto, x=x_axis, y="p99", hover_data=hover_fields, color=color_by, markers=True
        )
    fig.update_layout(title=title_text, title_x=0.5)
    fig.update_layout(plot_bgcolor="white")
    fig.update_xaxes(
        mirror=False, ticks="outside", showline=True, linecolor="black", gridcolor="lightgrey"
    )
    fig.update_yaxes(
        mirror=False, ticks="outside", showline=True, linecolor="black", gridcolor="lightgrey"
    )
    return fig
