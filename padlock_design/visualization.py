import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import io

def plot_structure(sequence, structure, title="Padlock Structure", regions=None):
    """
    Draws a secondary structure diagram using matplotlib and networkx.
    Uses a circular layout where backbone is the circle and pairs are chords.
    Highlights different regions of the backbone with colors.
    """
    sequence = sequence.upper()
    n = len(sequence)
    
    # Create graph and add nodes
    G = nx.Graph()
    base_colors = {'A': '#FF5252', 'T': '#448AFF', 'U': '#448AFF', 'G': '#69F0AE', 'C': '#FFAB40'}
    node_colors = [base_colors.get(base, '#E0E0E0') for base in sequence]
    for i, base in enumerate(sequence):
        G.add_node(i, label=base)
    
    # Add backbone edges with a default region
    for i in range(n - 1):
        G.add_edge(i, i + 1, type='backbone', region='Backbone') # Default region

    # Parse dot-bracket structure for pairing
    stack = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                start = stack.pop()
                G.add_edge(start, i, type='pair')

    # Assign specific regions to backbone edges if provided
    if regions:
        for name, (start, end) in regions.items():
            # Ensure indices are within bounds
            if start is None or end is None: continue
            start = max(0, start)
            end = min(n - 1, end)
            for i in range(start, end):
                if G.has_edge(i, i + 1):
                    G.edges[i, i + 1]['region'] = name

    # Layout
    pos = nx.circular_layout(G)
    
    # Draw
    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    
    # Define region colors
    region_colors = {
        "5' miRNA Arm": 'blue',
        "3' miRNA Arm": 'blue',
        'Cas12a binding site': 'green',
        'Spacer': 'purple',
        'Short spacer': 'orange',
        'Backbone': 'gray'
    }
    
    # Separate and draw backbone edges by region for coloring
    edges_by_region = {name: [] for name in region_colors}
    backbone_edges = [(u, v, d) for u, v, d in G.edges(data=True) if d.get('type') == 'backbone']
    for u, v, data in backbone_edges:
        region = data.get('region', 'Backbone')
        edges_by_region[region].append((u, v))

    for region, edge_list in edges_by_region.items():
        if edge_list:
            nx.draw_networkx_edges(G, pos, edgelist=edge_list, width=3.0, 
                                   edge_color=region_colors.get(region, 'gray'), alpha=0.8)

    # Draw pairing edges (pink dotted lines)
    pair_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('type') == 'pair']
    nx.draw_networkx_edges(G, pos, edgelist=pair_edges, width=1.5, edge_color='#FF4081', style='dashed')
    
    # Draw all nodes with a standard black border
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=500, edgecolors='black', linewidths=1.0)
    
    # Draw base labels inside nodes
    labels = {i: base for i, base in enumerate(sequence)}
    nx.draw_networkx_labels(G, pos, labels, font_size=10, font_weight='bold', font_color='black')
    
    # Annotate 5'-P at the start of the padlock
    if 0 in pos:
        x, y = pos[0]
        plt.text(x * 1.15, y, "5'-P", ha='left', va='center', fontsize=10, color='black', weight='bold')

    # Annotate indices (e.g., every 10th) to reduce clutter
    for i, (x, y) in pos.items():
        if i == 0: continue # Skip '1' to avoid overlap with 5'-P
        if (i + 1) % 10 == 0:
            dist = 1.15
            lbl_x = x * dist
            lbl_y = y * dist
            plt.text(lbl_x, lbl_y, str(i + 1), ha='center', va='center', fontsize=8, color='gray')

    # Add a legend for regions, positioned outside the plot
    handles = [plt.Line2D([0], [0], color=color, lw=4, label=name) 
               for name, color in region_colors.items() if name != 'Backbone' and edges_by_region.get(name)]
    if handles:
        ax.legend(handles=handles, title="Regions", bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)

    plt.title(title, fontsize=15, fontweight='bold', pad=20)
    plt.axis('off')
    plt.tight_layout(rect=[0, 0, 0.85, 1]) # Adjust layout to make room for legend
    
    # Save to buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight', dpi=150)
    plt.close()
    buf.seek(0)
    return buf
