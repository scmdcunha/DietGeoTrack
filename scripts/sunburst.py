#!/usr/bin/env python3
import pandas as pd
import plotly.express as px
from pathlib import Path

csv_file = Path("results/blast/arthropoda_top3_taxonomy.csv")

df = pd.read_csv(csv_file)

# Counts the number of hits per combination of Order -> Family -> Genus
count_df = (
    df
    .groupby(["Order", "Family", "Genus"])
    .size()
    .reset_index(name="Count")
)

# Creates the sunburst:
#  - path defines hierarchy
#  - values defines the size of slices
#  - color can be Order to give distinct colors to orders
fig = px.sunburst(
    count_df,
    path=["Order", "Family", "Genus"],
    values="Count",
    color="Order",
    title="Sunburst: Orders → Families → Genera in Bat Diet"
)

fig.update_layout(margin=dict(t=40, l=0, r=0, b=0))

# Saves as interactive HTML
output_html = "results/blast/taxonomy_sunburst.html"
fig.write_html(output_html)

print(f"Sunburst saved to {output_html}")
