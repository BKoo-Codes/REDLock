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
    
    # Add backbone edges
    for i in range(n - 1):
        G.add_edge(i, i + 1, type='backbone')
    # Padlock is circular DNA, so connect ends? 
    # Usually padlock probes are linear until ligated. Let's assume linear display for now,
    # or circular if the user considers it ligated. But input is usually linear probe.
    # Let's draw it as a circle to show the "potential" structure clearly.
    G.add_edge(n - 1, 0, type='backbone') # Connect ends for visual circularity

    # Parse dot-bracket structure to add pairing edges
    stack = []
    pairs = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                start = stack.pop()
                pairs.append((start, i))
                G.add_edge(start, i, type='pair')

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
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=500, edgecolors='black', linewidths=1.0)
    
    # Draw labels
    labels = {i: base for i, base in enumerate(sequence)}
    nx.draw_networkx_labels(G, pos, labels, font_size=10, font_weight='bold', font_color='black')
    
    # Annotate indices (every 5 or 10)
    # Using a slightly larger radius for index labels
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
