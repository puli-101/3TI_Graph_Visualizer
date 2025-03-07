import networkx as nx
import matplotlib.pyplot as plt
import os 

"""
Graph display and image/graph serialization functions
"""

#Saves displayed graph into png image
def save_graph(n,m,k,q):
    #n,m,k dimensions of graph
    #q order of field

    # Define the directory and base file name
    directory = '/home/david/Desktop/3TIP/'
    base_filename = f'T_n{n}_m{m}_k{k}_q{q}'
    extension = '.png'

    # Ensure the directory exists
    os.makedirs(directory, exist_ok=True)

    # Initialize a counter and check if the file exists
    n = 0
    while os.path.exists(os.path.join(directory, f"{base_filename}_{n}{extension}")):
        n += 1

    # Create the new file name
    output_path = os.path.join(directory, f"{base_filename}_{n}{extension}")
    #dpi_value = 300

    plt.savefig(output_path, bbox_inches='tight') #dpi=dpi_value


#Displays sage graph by translating it into NX
def graph_display(G,n,m,k,q):
    #G sage graph
    #n,m,k dimensions
    #q field

    G_vis = nx.Graph()
    # Add nodes and edges
    G_vis.add_nodes_from(G.vertices())
    G_vis.add_edges_from([(u,v) for u,v,_ in G.edges()])

    # Draw the graph
    nx.draw(G_vis, with_labels=False, node_color='white', edgecolors='black')
    save_graph(n,m,k,q)
    plt.show()