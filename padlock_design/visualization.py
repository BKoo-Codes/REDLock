import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
import networkx as nx
import numpy as np
import io

# ===== Dark / neon fluorescence palette (matches the app UI) =====
BASE_COLORS = {
    'A': '#ff6f6f',   # coral
    'T': '#22d3ee',   # cyan
    'U': '#22d3ee',   # cyan
    'G': '#3ef2b3',   # signal green
    'C': '#9b7bff',   # violet
}
COL_BACKBONE = '#7d8bc4'
COL_PAIR     = '#ec5fd0'   # magenta base-pair bonds
COL_TEXT     = '#eaf0ff'
COL_MUTED    = '#93a0c4'
COL_NODE_EDGE = '#0b0f20'
COL_LABEL    = '#06121a'   # dark text on bright neon nodes


def plot_structure(sequence, structure, title="Padlock Structure"):
    """
    Draws a secondary structure diagram on a transparent dark canvas so it blends
    into the app's dark theme. Uses a circular layout with a visible 5'/3' gap,
    neon base colouring, glowing base-pair bonds, and a compact base legend.
    """
    sequence = sequence.upper()
    n = len(sequence)

    # Build graph
    G = nx.Graph()
    node_colors = [BASE_COLORS.get(base, '#c4d2db') for base in sequence]
    for i, base in enumerate(sequence):
        G.add_node(i, label=base)

    # Dummy node creates a visual gap between start and end in the circular layout
    G.add_node('gap_node')

    # Backbone edges
    for i in range(n - 1):
        G.add_edge(i, i + 1, type='backbone')

    # Parse dot-bracket pairing
    stack = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                start = stack.pop()
                G.add_edge(start, i, type='pair')

    pos = nx.circular_layout(G)

    # Canvas (transparent so the dark app card shows through)
    fig = plt.figure(figsize=(8, 8))
    ax = plt.gca()
    fig.patch.set_alpha(0.0)
    ax.set_facecolor('none')

    backbone_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('type') == 'backbone']
    pair_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('type') == 'pair']

    # Base-pair bonds: soft glow underneath + crisp dashed line on top
    if pair_edges:
        nx.draw_networkx_edges(G, pos, edgelist=pair_edges, width=6.0,
                               edge_color=COL_PAIR, alpha=0.18)
        nx.draw_networkx_edges(G, pos, edgelist=pair_edges, width=1.8,
                               edge_color=COL_PAIR, style='dashed', alpha=0.95)

    # Backbone
    nx.draw_networkx_edges(G, pos, edgelist=backbone_edges, width=3.0,
                           edge_color=COL_BACKBONE, alpha=0.9)

    # Nodes: glow halo then solid neon base
    actual_nodes = list(range(n))
    node_col = [node_colors[i] for i in actual_nodes]
    nx.draw_networkx_nodes(G, pos, nodelist=actual_nodes, node_color=node_col,
                           node_size=900, alpha=0.20)
    nx.draw_networkx_nodes(G, pos, nodelist=actual_nodes, node_color=node_col,
                           node_size=480, edgecolors=COL_NODE_EDGE, linewidths=1.2)

    # Base letters (dark text reads well on bright neon fills)
    labels = {i: base for i, base in enumerate(sequence)}
    nx.draw_networkx_labels(G, pos, labels, font_size=10, font_weight='bold',
                            font_color=COL_LABEL)

    # Index annotations (1, milestones every 10, and n)
    dist = 1.12
    fontsize = 12
    for i in range(n):
        (x_node, y_node) = pos[i]
        lbl_x, lbl_y = x_node * dist, y_node * dist
        if i == 0:
            plt.text(lbl_x, lbl_y, "1", ha='center', va='center',
                     fontsize=fontsize, color=COL_MUTED, weight='bold')
        elif i == n - 1:
            plt.text(lbl_x, lbl_y, str(n), ha='center', va='center',
                     fontsize=fontsize, color=COL_MUTED, weight='bold')
        elif (i + 1) % 10 == 0:
            plt.text(lbl_x, lbl_y, str(i + 1), ha='center', va='center',
                     fontsize=fontsize, color=COL_MUTED)

    # Title with a subtle dark glow so it stays readable
    t = plt.title(title, fontsize=24, fontweight='bold', color=COL_TEXT, pad=20)
    t.set_path_effects([pe.withStroke(linewidth=3, foreground='#050810')])

    # Compact base legend
    handles = [Line2D([0], [0], marker='o', linestyle='None',
                      markerfacecolor=BASE_COLORS[b], markeredgecolor=COL_NODE_EDGE,
                      markersize=11, label=b) for b in ['A', 'T', 'G', 'C']]
    ax.legend(handles=handles, loc='upper right', frameon=False, fontsize=11,
              labelcolor=COL_TEXT, handletextpad=0.3, ncol=4, columnspacing=0.9,
              bbox_to_anchor=(1.0, 1.03))

    plt.axis('off')
    plt.tight_layout()

    # Save (transparent PNG blends into the dark app)
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight', dpi=150, transparent=True)
    plt.close()
    buf.seek(0)
    return buf
