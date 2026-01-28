import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import io

def plot_structure(sequence, structure, title="Padlock Structure", regions=None):
    """
    Draws a secondary structure diagram using matplotlib and networkx.
    Uses a circular layout with a visible gap between start and end.
    Highlights different regions of the backbone with colors and a legend.
    """
    sequence = sequence.upper()
    n = len(sequence)
    
    # Create graph and add nodes
    G = nx.Graph()
    base_colors = {'A': '#FF5252', 'T': '#448AFF', 'U': '#448AFF', 'G': '#69F0AE', 'C': '#FFAB40'}
    node_colors = [base_colors.get(base, '#E0E0E0') for base in sequence]
    for i, base in enumerate(sequence):
        G.add_node(i, label=base)
    
    # --- Add a dummy node to create a visual gap in the circular layout ---
    G.add_node('gap_node')
    
    # Add backbone edges with a default region
    for i in range(n - 1):
        G.add_edge(i, i + 1, type='backbone', region='Backbone')

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
            if start is None or end is None: continue
            start = max(0, start)
            end = min(n, end) # Use n for end range
            for i in range(start, end):
                if i < n -1 and G.has_edge(i, i + 1):
                    G.edges[i, i + 1]['region'] = name

    # --- Layout: Use circular_layout with the dummy node to create a gap ---
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
    
    # Separate and draw backbone edges by region
    edges_by_region = {name: [] for name in region_colors}
    backbone_edges = [(u, v, d) for u, v, d in G.edges(data=True) if d.get('type') == 'backbone']
    for u, v, data in backbone_edges:
        region = data.get('region', 'Backbone')
        edges_by_region[region].append((u, v))

    for region, edge_list in edges_by_region.items():
        if edge_list:
            nx.draw_networkx_edges(G, pos, edgelist=edge_list, width=3.5, 
                                   edge_color=region_colors.get(region, 'gray'), alpha=0.9)

    # Draw pairing edges
    pair_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('type') == 'pair']
    nx.draw_networkx_edges(G, pos, edgelist=pair_edges, width=1.5, edge_color='#FF4081', style='dashed')
    
    # --- Draw only the actual nodes, not the dummy gap_node ---
    actual_nodes = list(range(n))
    nx.draw_networkx_nodes(G, pos, nodelist=actual_nodes, node_color=[node_colors[i] for i in actual_nodes], 
                           node_size=500, edgecolors='black', linewidths=1.0)
    
    # Draw base labels inside nodes
    labels = {i: base for i, base in enumerate(sequence)}
    nx.draw_networkx_labels(G, pos, labels, font_size=10, font_weight='bold', font_color='black')
    
    # --- Annotation for Indices and Termini ---
    for i in range(n):
        (x, y) = pos[i]
        is_first = (i == 0)
        is_last = (i == n - 1)
        is_milestone = (i + 1) % 10 == 0 and not is_first and not is_last

        if is_first or is_last or is_milestone:
            dist = 1.18
            lbl_x, lbl_y = x * dist, y * dist
            
            if is_first:
                plt.text(x * 1.25, y, "5'-P", ha='left', va='center', fontsize=10, color='black', weight='bold')
                plt.text(x * 1.1, y, "1", ha='right', va='center', fontsize=8, color='gray', weight='bold')
            elif is_last:
                plt.text(lbl_x, lbl_y, str(n), ha='center', va='center', fontsize=8, color='gray', weight='bold')
            else: # Milestones
                plt.text(lbl_x, lbl_y, str(i + 1), ha='center', va='center', fontsize=8, color='gray')

    # Add a legend for regions
    handles = [plt.Line2D([0], [0], color=color, lw=4, label=name) 
               for name, color in region_colors.items() if name != 'Backbone' and edges_by_region.get(name)]
    if handles:
        ax.legend(handles=handles, title="Regions", bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)

    plt.title(title, fontsize=15, fontweight='bold', pad=20)
    plt.axis('off')
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    
    # Save to buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight', dpi=150)
    plt.close()
    buf.seek(0)
    return buf
