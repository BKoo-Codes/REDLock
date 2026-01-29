import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import io

def plot_structure(sequence, structure, title="Padlock Structure"):
    """
    Draws a secondary structure diagram using matplotlib and networkx.
    Uses a circular layout with a visible gap between start and end.
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
    
    # Add backbone edges
    for i in range(n - 1):
        G.add_edge(i, i + 1, type='backbone')

    # Parse dot-bracket structure for pairing
    stack = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                start = stack.pop()
                G.add_edge(start, i, type='pair')

    # --- Layout: Use circular_layout with the dummy node to create a gap ---
    pos = nx.circular_layout(G)
    
    # Draw
    plt.figure(figsize=(8, 8))
    ax = plt.gca()
    
    # --- Draw backbone edges in a single color ---
    backbone_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('type') == 'backbone']
    nx.draw_networkx_edges(G, pos, edgelist=backbone_edges, width=3.0, 
                           edge_color='gray', alpha=0.9)

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
    
    # --- Annotation for Indices ---
    dist = 1.12  # Distance for number labels
    fontsize = 9 # Font size for number labels
    
    for i in range(n):
        (x_node, y_node) = pos[i] # Position of the node
        is_first = (i == 0)
        is_last = (i == n - 1)
        is_milestone = (i + 1) % 10 == 0

        # Calculate position for number labels
        lbl_x, lbl_y = x_node * dist, y_node * dist
        
        # Draw number labels
        if is_first:
            # Draw '1'
            plt.text(lbl_x, lbl_y, "1", ha='center', va='center', fontsize=fontsize, color='gray', weight='bold')
        elif is_last:
            # Draw 'n'
            plt.text(lbl_x, lbl_y, str(n), ha='center', va='center', fontsize=fontsize, color='gray', weight='bold')
        elif is_milestone:
            plt.text(lbl_x, lbl_y, str(i + 1), ha='center', va='center', fontsize=fontsize, color='gray')

    plt.title(title, fontsize=15, fontweight='bold', pad=20)
    plt.axis('off')
    plt.tight_layout()
    
    # Save to buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight', dpi=150)
    plt.close()
    buf.seek(0)
    return buf
