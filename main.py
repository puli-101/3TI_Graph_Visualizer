import time
import argparse
from tensor import *
from tools import *

def argparser():
    #Parses the values of (n,m,k,q,labeled) as described
    parser = argparse.ArgumentParser(
        description="Build graph from a random 3-tensor",
        epilog="Example usage: sage script.sage -n 4 -m 3 -k 5 -q 7 -labeled -verbose --isolated_nodes"
    )
    parser.add_argument("-n", type=int, default=4, help="Dimension n for the first vector space")
    parser.add_argument("-m", type=int, default=4, help="Dimension m for the second vector space")
    parser.add_argument("-k", type=int, default=4, help="Dimension k for the third vector space")
    parser.add_argument("-q", type=int, default=5, help="Prime field size")
    parser.add_argument("--labeled", action="store_true", help="Show graph with vertex labels")
    parser.add_argument("--verbose", action="store_true", help="Show extra info on terminal")
    parser.add_argument("--isolated_nodes", action="store_true", help="Displays nodes of degree zero on the final graph")
    parser.add_argument("--degree_ubound", type=int, default=100, help="Filters all nodes of degree greater or equal than specified")
    parser.add_argument("--degree_lbound", type=int, default=0, help="Filters all nodes of degree lower or equal than specified")

    return parser.parse_args()

if __name__ == "__main__":
    args = argparser()
    
    #Dimensions of U, V, W.
    n = args.n
    m = args.m
    k = args.k
    #Field size
    q = args.q
    #Field
    F = GF(q)  

    #Show graph labels?
    labeled = args.labeled
    #Show extra info?
    verbose = args.verbose
    #Display degree zero nodes?
    deg_0 = args.isolated_nodes 

    #check passed parameters
    print(n,m,k,q,labeled, verbose)
    
    #Create a random 3-tensor T of dimensions n x m x k over F.
    T = [[[F.random_element() for _ in range(k)] for _ in range(m)] for _ in range(n)]
    start = time.time()
    G = tensor_to_graph(T, n, m, k, F)
    print(f"Computation time: {time.time() - start}")
    print("Tensor T:")
    for i in range(n):
        print(T[i])
    
    if verbose:
        print("\nGraph vertices: ")
        print(G.vertices())
        print("\nGraph edges: ")
        print(G.edges())

    #Identify isolated vertices
    if not(deg_0):
        print("Removing isolated vertices")
        isolated_vertices = [v for v in G.vertices() if G.degree(v) == 0]
        G.delete_vertices(isolated_vertices)

    #Display graph
    graph_display(G,n,m,k,q)
    