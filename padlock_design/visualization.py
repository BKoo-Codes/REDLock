import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import io

def plot_structure(sequence, structure, title="Padlock Structure"):
    """
    Draws a secondary structure diagram using matplotlib and networkx.
    Uses a circular layout where backbone is the circle and pairs are chords.
    """
    sequence = sequence.upper()
    n = len(sequence)
    
    # Create graph
    G = nx.Graph()
    
    # Add nodes (bases)
    base_colors = {'A': '#FF5252', 'T': '#448AFF', 'U': '#448AFF', 'G': '#69F0AE', 'C': '#FFAB40'}
    node_colors = [base_colors.get(base, '#E0E0E0') for base in sequence]
    
    for i, base in enumerate(sequence):
        G.add_node(i, label=base)
    
    # Add backbone edges, but not between the last and first node to create a gap.
    for i in range(n - 1):
        G.add_edge(i, i + 1, type='backbone')

    # Parse dot-bracket structure to add pairing edges
    stack = []
    pairs = []
    binding_nodes_set = set()
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                start = stack.pop()
                pairs.append((start, i))
                G.add_edge(start, i, type='pair')
                binding_nodes_set.add(start)
                binding_nodes_set.add(i)

    # Layout: Circular
    pos = nx.circular_layout(G)
    
    # Draw
    plt.figure(figsize=(8, 8))
    ax = plt.gca()
    
    # Draw backbone edges
    backbone_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('type') == 'backbone']
    nx.draw_networkx_edges(G, pos, edgelist=backbone_edges, width=2.0, edge_color='gray', alpha=0.5)
    
    # Draw pair edges (chords)
    pair_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('type') == 'pair']
    nx.draw_networkx_edges(G, pos, edgelist=pair_edges, width=1.5, edge_color='#FF4081', style='dashed', alpha=0.8)
    
    # Get all nodes and separate them
    all_nodes = list(G.nodes())
    binding_nodes = list(binding_nodes_set)
    non_binding_nodes = [node for node in all_nodes if node not in binding_nodes_set]

    # Draw non-binding nodes
    nx.draw_networkx_nodes(G, pos, nodelist=non_binding_nodes, 
                           node_color=[node_colors[i] for i in non_binding_nodes],
                           node_size=500, edgecolors='black', linewidths=1.0)
    
    # Draw binding nodes with a dashed border
    if binding_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=binding_nodes, 
                               node_color=[node_colors[i] for i in binding_nodes],
                               node_size=500, edgecolors='#FF4081', linewidths=2.0, linestyle='dashed')
    
    # Draw labels
    labels = {i: base for i, base in enumerate(sequence)}
    nx.draw_networkx_labels(G, pos, labels, font_size=10, font_weight='bold', font_color='black')
    
    # Annotate 5'-P at the start of the padlock
    if 0 in pos:
        x, y = pos[0]
        offset = 0.15
        plt.text(x + offset, y, "5'-P", ha='left', va='center', fontsize=10, color='black', weight='bold')

    # Annotate indices (every 5 or 10)
    for i, (x, y) in pos.items():
        if (i + 1) % 5 == 0 or i == 0:
            angle = np.arctan2(y, x)
            dist = 1.15 # Distance factor
            lbl_x = x * dist
            lbl_y = y * dist
            plt.text(lbl_x, lbl_y, str(i + 1), ha='center', va='center', fontsize=8, color='gray')

    plt.title(title, fontsize=15, fontweight='bold', pad=20)
    plt.axis('off')
    
    # Save to buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight', dpi=150)
    plt.close()
    buf.seek(0)
    return buf
