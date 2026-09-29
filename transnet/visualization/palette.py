"""One colour scheme for every TransNet figure, validated rather than eyeballed.

Colour does exactly two jobs in these figures, and they never share a hue:

**Direction** -- up/down, activating/inhibiting -- is *polarity*, so it takes a
diverging pair: red for increase and activation, blue for decrease and
inhibition, grey for "measured, no change". Activation shares red with increase
because an activator's effect on a reaction is an increase.

**Identity** -- which layer, which regulation axis, which concordance class --
is *categorical*, so it takes hues that are neither red nor blue. The previous
scheme (matplotlib's tab10) coloured the Reactions layer red and the
Transcriptome blue, so a red ring and a red fill on the same node meant
different things.

Validation
----------
Checked with the dataviz palette validator (OKLab Delta E x100 under simulated
protanopia and deuteranopia, Machado 2009; normal-vision floor 15; lightness
band and chroma floor) on the light chart surface ``#fcfcfb``:

* Layer set (Signaling, Transcriptome, Proteome, Reactions, Metabolome),
  **all pairs**: passes every hard gate. Worst CVD pair Reactions/Metabolome
  Delta E 6.1 (deutan) sits in the 6-8 band, which is legal only with
  secondary encoding -- every figure provides it: Reactions are squares, every
  plane and bar is labelled by name.
* Direction pair (blue/red), all pairs: CVD Delta E 21.6, normal 32.3.
* Yellow, aqua and magenta are below 3:1 contrast on the surface, so they are
  never used for text and always sit beside a visible label.

Text never takes a series colour; it uses the ink tokens below.
"""

__all__ = [
    "LAYER_COLORS", "UP", "DOWN", "SIGN_COLORS", "UNCHANGED", "UNMEASURED",
    "UNDIRECTED", "AXIS_COLORS", "CONCORDANCE_COLORS", "INK", "SURFACE",
    "GRID", "BASELINE", "MUTED", "SECONDARY", "FLOW", "style_axes",
]

# -- identity: layers ------------------------------------------------------------
LAYER_COLORS = {
    "Signaling": "#008300",      # green
    "Transcriptome": "#4a3aa7",  # violet
    "Proteome": "#eda100",       # yellow
    "Reactions": "#1baf7a",      # aqua
    "Metabolome": "#e87ba4",     # magenta
    "Pathways": "#898781",
}

# -- direction -------------------------------------------------------------------
UP = "#e34948"
DOWN = "#2a78d6"
UNCHANGED = "#c3c2b7"      # measured, no significant change: the neutral midpoint
UNMEASURED = "#fcfcfb"     # surface: a hollow node
UNDIRECTED = "#52514e"     # changed, direction unknown (omnibus tests)
SIGN_COLORS = {1: UP, -1: DOWN, 0: UNCHANGED}

# -- identity derived from layers -------------------------------------------------
#: A regulation axis is coloured by the layer that carries it; the two axes
#: meet at the reaction, so "both" takes the Reactions colour.
AXIS_COLORS = {
    "enzyme": LAYER_COLORS["Proteome"],
    "metabolite": LAYER_COLORS["Metabolome"],
    "both": LAYER_COLORS["Reactions"],
}

CONCORDANCE_COLORS = {
    "concordant": LAYER_COLORS["Reactions"],
    "protein_only": LAYER_COLORS["Proteome"],
    "transcript_only": LAYER_COLORS["Transcriptome"],
    "discordant": LAYER_COLORS["Metabolome"],
    "unchanged": UNCHANGED,
}

# -- ink and chrome -----------------------------------------------------------------
INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
SURFACE = "#fcfcfb"
FLOW = "#898781"


def style_axes(axes, grid_axis: str = "x"):
    """Recessive chrome: hairline solid grid on one axis, no top/right spines."""
    axes.set_facecolor(SURFACE)
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_color(BASELINE)
        axes.spines[side].set_linewidth(0.8)
    axes.tick_params(colors=SECONDARY, labelcolor=SECONDARY, length=3, width=0.6)
    if grid_axis:
        axes.grid(axis=grid_axis, color=GRID, linewidth=0.6, linestyle="-")
        axes.set_axisbelow(True)
    axes.title.set_color(INK)
    axes.xaxis.label.set_color(SECONDARY)
    axes.yaxis.label.set_color(SECONDARY)
    return axes
