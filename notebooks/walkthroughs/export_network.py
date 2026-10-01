# %% [markdown]
# # Saving, exporting and sharing a network
#
# A trans-omic network is usually needed outside Python: in Cytoscape for a
# figure, in a 3-D viewer for a talk, as files for a collaborator, or as a
# matrix for another tool. This walkthrough writes the responsive network of
# the bundled example in every format TransNet supports and explains who each
# one is for.
#
# | Format | Function | For |
# |---|---|---|
# | Two CSV files | `write_network` / `read_network` | archiving a network and reloading it in TransNet |
# | Two CSV files, from the builder | `Transnet.save_network` / `Transnet.load_network` | networks built with the `Transnet` object API |
# | Adjacency matrix | `Transnet.generate_adjacency_matrix` | matrix-based tools (R, MATLAB, spreadsheets) |
# | Cytoscape JSON | `to_cytoscape_json` | Cytoscape Desktop and the cytoscape.js web library |
# | Arena3D Web JSON | `to_arena3d` | interactive 3-D multilayer viewing at arena3d.org |
# | transomics2cytoscape bundle | `to_transomics2cytoscape` | stacked 3-D layer figures in Cytoscape, driven from R |
# | Self-contained HTML | `plot_transomic_network_interactive` | sending an interactive figure to someone without Python |
# | Static figure (PNG, PDF, SVG) | any `plot_*` function, then `savefig` | papers and slides |
#
# Everything is written to a folder called `exports/` in the working
# directory.

# %%
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from transnet import (
    Transnet,
    available_edge_types,
    load_example_network,
    load_example_omics,
    map_omics_to_network,
    responsive_subnetwork,
)
from transnet.io import (
    read_network,
    to_arena3d,
    to_cytoscape_json,
    to_transomics2cytoscape,
    write_network,
)
from transnet.visualization import (
    EDGE_STYLES,
    LAYER_COLORS,
    plot_transomic_network,
    plot_transomic_network_interactive,
)

graph = load_example_network()
map_omics_to_network(
    graph, load_example_omics(), id_column="id", log2fc_column="log2FC", qvalue_column="padj",
)
responsive = responsive_subnetwork(graph)

out = Path("exports")
out.mkdir(exist_ok=True)
print(f"responsive network: {responsive.number_of_nodes()} molecules, "
      f"{responsive.number_of_edges()} relationships")

# %% [markdown]
# ## CSV: the TransNet format
#
# `write_network` writes two files: `nodes.csv` (one row per molecule, with its
# layer, name and measured values) and `interactions.csv` (one row per edge,
# with its type, sign and source database). These are plain tables that open
# in any spreadsheet, and `read_network` loads them back with every attribute
# the analyses need. The organism networks in `data/<organism>/latest/` are
# stored the same way.

# %%
write_network(responsive, str(out / "csv"))
print(", ".join(sorted(p.name for p in (out / "csv").iterdir())))

reloaded = read_network(str(out / "csv" / "interactions.csv"),
                        nodes_file=str(out / "csv" / "nodes.csv"))
print(f"reloaded {reloaded.number_of_nodes()} molecules and "
      f"{reloaded.number_of_edges()} relationships; edge types preserved: "
      f"{available_edge_types(reloaded) == available_edge_types(responsive)}")
pd.read_csv(out / "csv" / "interactions.csv").head(4)

# %% [markdown]
# ## CSV from the builder object
#
# A network assembled with the `Transnet` object API (see the
# *trans-omic network* walkthrough) has its own `save_network` and
# `load_network` methods, which write the same two files from the layer
# objects. `save_adjacency=True` also writes the adjacency matrix. Here a small
# builder network is created from one reaction so the round trip is quick.

# %%
from transnet import Metabolite, Protein, Reaction
from transnet.biology.layers import Metabolome, Proteome, Reactions

network = Transnet(name="hexokinase")
hexokinase = Protein(uniprot_id="P52789", name="Hexokinase-2", gene=["Hk2"],
                     ec_number=["2.7.1.1"])
network.proteome = Proteome()
network.proteome.proteins = [hexokinase]
network.metabolome = Metabolome()
network.metabolome.metabolites = [
    Metabolite(kegg_compound_id="C00031", kegg_name="D-Glucose"),
    Metabolite(kegg_compound_id="C00092", kegg_name="D-Glucose 6-phosphate"),
]
network.reactions = Reactions()
network.reactions.reactions = [Reaction(
    id="R00299", name="hexokinase", enzyme=["2.7.1.1"],
    substrates=["C00031"], products=["C00092"], reversible=False,
)]
network.generate_graph()

network.save_network(str(out / "builder"), save_adjacency=True)
print(", ".join(sorted(p.name for p in (out / "builder").iterdir())))
again = Transnet.load_network(str(out / "builder"))
print(f"loaded back: {len(again.proteome.proteins)} protein, "
      f"{len(again.metabolome.metabolites)} metabolites, "
      f"{len(again.reactions.reactions)} reaction")

# %% [markdown]
# ## Adjacency matrix
#
# `generate_adjacency_matrix` returns the network as a square table: one row and
# one column per molecule, with a non-zero entry where an edge joins them.
# `symmetric=False` keeps the direction, so entry (i, j) is the edge from i to
# j. A dense matrix grows with the square of the number of nodes, so use it only
# for small networks such as a responsive subnetwork.

# %%
matrix = network.generate_adjacency_matrix(symmetric=False)
matrix

# %% [markdown]
# ## Cytoscape JSON
#
# `to_cytoscape_json` writes the format that Cytoscape Desktop reads through
# *File > Import > Network from File* and that the cytoscape.js JavaScript
# library reads directly. Every node keeps its layer, name and fold change, and
# every edge its type and sign, so Cytoscape's style editor can colour and shape
# them from those columns.

# %%
to_cytoscape_json(responsive, out / "network_cytoscape.json")
elements = json.loads((out / "network_cytoscape.json").read_text())["elements"]
print(f"{len(elements['nodes'])} nodes and {len(elements['edges'])} edges")
elements["nodes"][0]["data"]

# %% [markdown]
# ## Arena3D Web
#
# [Arena3D Web](https://arena3d.pavlopouloslab.info/) draws each layer as a
# plane in 3-D, which makes the edges that cross layers easy to see. Load the
# file through Arena3D Web's session upload. `max_nodes_per_layer` keeps the file small
# enough for the browser by dropping the least connected nodes.

# %%
to_arena3d(responsive, out / "network_arena3d.json", max_nodes_per_layer=500)
session = json.loads((out / "network_arena3d.json").read_text())
print("top-level keys:", ", ".join(session))

# %% [markdown]
# ## transomics2cytoscape
#
# [transomics2cytoscape](https://bioconductor.org/packages/transomics2cytoscape/)
# is a Bioconductor package that draws trans-omic networks in Cytoscape as
# stacked layers, the style used in the Kuroda laboratory's papers.
# `to_transomics2cytoscape` writes one node table and one edge table per layer,
# plus a table of the edges between layers, and a `README.md` with the R code
# that loads them. It returns the files as a zip archive, or writes them to a
# folder with `as_zip=False`.

# %%
to_transomics2cytoscape(responsive, zip_path=out / "network_transomics2cytoscape.zip")
to_transomics2cytoscape(responsive, output_dir=out / "transomics2cytoscape", as_zip=False)
sorted(str(p.relative_to(out / "transomics2cytoscape"))
       for p in (out / "transomics2cytoscape").rglob("*") if p.is_file())

# %% [markdown]
# ## A self-contained interactive page
#
# `plot_transomic_network_interactive` returns a plotly figure, and
# `write_html` saves it as a single HTML file that opens in any browser with no
# Python and no internet connection. This is usually the best way to send a
# network to a collaborator. `layout="layered"` stacks the layers in rows, as
# in the static figures.

# %%
figure = plot_transomic_network_interactive(responsive, layout="layered",
                                            title="Responsive trans-omic network")
figure.write_html(out / "network.html", include_plotlyjs=True)
print(f"network.html: {(out / 'network.html').stat().st_size / 1e6:.1f} MB")

# %% [markdown]
# The same figure inline. Hover over a node to see what it is and how it
# changed, drag to zoom, and click a legend entry to hide a layer or an edge
# type.

# %%
import plotly.io as pio

pio.renderers.default = "notebook_connected"
figure

# %% [markdown]
# ## Static figures
#
# Every `plot_*` function returns a matplotlib figure, so `savefig` writes it in
# any format matplotlib supports: PNG for slides, and PDF or SVG for papers,
# where the text stays editable.

# %%
static = plot_transomic_network(responsive, title="Responsive trans-omic network")
for suffix in ("png", "pdf", "svg"):
    static.savefig(out / f"network.{suffix}", dpi=200, bbox_inches="tight")
plt.show()

# %% [markdown]
# The colours and line styles used in every TransNet figure are available as
# `LAYER_COLORS` and `EDGE_STYLES`, so a figure made with another tool can
# match them.

# %%
pd.Series(LAYER_COLORS, name="colour").to_frame()

# %%
pd.Series(EDGE_STYLES, name="line style").to_frame()

# %% [markdown]
# ## Everything written

# %%
sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_file())
