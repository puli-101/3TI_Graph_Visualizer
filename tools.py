import networkx as nx
import matplotlib.pyplot as plt
import os 

"""
Graph display and image/graph serialization functions
"""


#Finds all unique simple cycles of length c in a graph G
#used in graph display to color these ndoes differently
def find_cycles_of_length_c(G, c):
    """
    G: A SageMath graph object
    c: The desired cycle length

    return: A list of cycles, where each cycle is represented as a tuple of vertices

    A cycle is considered unique up to cyclic permutation and reversal.
    To avoid duplicates, we only add a cycle if its starting vertex is the smallest 
    in the cycle
    """
    cycles_set = set()
    
    def dfs(start, current, path, visited):
        if len(path) > c:
            return

        if len(path) == c:
            #check if there's an edge from the current vertex back to the start
            if start in G.neighbors(current):
                cycle = path[:]  # A candidate cycle
                #to avoid over populating cycles_set only add cycle if start node has minimal label
                if start == min(cycle):
                    #At this point one could create a canonical representation:
                    #Compare the ordering of the cycle (after the start) with its reverse.
                    #forward = cycle[1:]
                    #backward = list(reversed(cycle[1:]))
                    #if backward < forward:
                    #    canonical_cycle = (cycle[0],) + tuple(backward)
                    #else:
                    #    canonical_cycle = tuple(cycle)
                    #cycles_set.add(canonical_cycle)
                    #however it does not really matter in terms of computation
                    cycles_set.update(cycle)
            return
        
        for neighbor in G.neighbors(current):
            if neighbor not in visited:
                dfs(start, neighbor, path + [neighbor], visited | {neighbor})
    
    #DFS starting from each node in graph
    for vertex in G.vertices():
        dfs(vertex, vertex, [vertex], {vertex})
        
    return list(cycles_set)


#Saves displayed graph into png image
def save_graph(n,m,k,q):
    #n,m,k dimensions of graph
    #q order of field

    #Define the directory and base file name
    directory = './graph_images/'
    base_filename = f'T_n{n}_m{m}_k{k}_q{q}'
    extension = '.png'

    #Ensure the directory exists
    os.makedirs(directory, exist_ok=True)

    #Initialize a counter and check if the file exists
    n = 0
    while os.path.exists(os.path.join(directory, f"{base_filename}_{n}{extension}")):
        n += 1

    #Create the new file name
    output_path = os.path.join(directory, f"{base_filename}_{n}{extension}")
    #dpi_value = 300

    plt.savefig(output_path, bbox_inches='tight') #dpi=dpi_value


#Displays sage graph by translating it into NX
def graph_display(G,n,m,k,q,cycle=None,labeled=False, save=False, loose=False, minimal=False):
    #G sage graph
    #n,m,k dimensions
    #q field

    #set fullscreen
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())

    #transform SageMath graph to NX graph
    G_vis = nx.Graph()

    #Add nodes and edges
    G_vis.add_nodes_from(G.vertices())
    G_vis.add_edges_from([(u,v) for u,v,_ in G.edges()])


    #if there is a specific type of cycle to compute
    if cycle != None and cycle > 2:
        u_bound = 3 if loose else cycle
        special_nodes = []
        # Define special nodes (e.g., nodes divisible by 20)
        for i in range(u_bound, cycle+1):
            print(f"Finding cycles of length {i}...")
            special_nodes += find_cycles_of_length_c(G,i)
            print(special_nodes)
        #Use a layout for consistent positioning
        pos = nx.spring_layout(G_vis)

        #Draw default nodes with a default color
        nx.draw_networkx_nodes(G_vis, pos, node_color='white', edgecolors='black')

        #Draw special nodes with a specific color
        nx.draw_networkx_nodes(G_vis, pos, nodelist=special_nodes, node_color='red', edgecolors='black')

        #Draw the edges and labels
        nx.draw_networkx_edges(G_vis, pos)
        #nx.draw_networkx_labels(G_vis, pos)
    else:
        #Otherwise just draw the graph
        nx.draw(G_vis, with_labels=labeled, node_color='white', edgecolors='black', node_size=8)

    if not(minimal):
        #save graph to file
        if save:
            save_graph(n,m,k,q)
        
        plt.show()

    