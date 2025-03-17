import time
import argparse
from tensor import *
from tools import *
from threading import Thread

def argparser():
    #Parses the values of (n,m,k,q,labeled) as described
    parser = argparse.ArgumentParser(
        description="Build graph from a random 3-tensor",
        epilog="Example usage: sage script.sage -n=4 -m=3 -k=5 -q=7 -c=3 --strict --labeled --verbose"
    )
    parser.add_argument("-n", type=int, default=4, help="Dimension n for the first vector space")
    parser.add_argument("-m", type=int, default=4, help="Dimension m for the second vector space")
    parser.add_argument("-k", type=int, default=4, help="Dimension k for the third vector space")
    parser.add_argument("-q", type=int, default=5, help="Prime field size")
    parser.add_argument("-c", type=int, default=None, help="Highlights all cycles of length 2 < c' <= c in the final graph")
    parser.add_argument("--strict", action="store_true", help="Only highlights cycle of length c")
    #parser.add_argument("--deg_ubound", type=int, default=1000, help="Filters all nodes of degree greater or equal than specified")
    #parser.add_argument("--deg_lbound", type=int, default=0, help="Filters all nodes of degree less or equal than specified")
    parser.add_argument("--labeled", action="store_true", help="Show graph with vertex labels")
    parser.add_argument("--verbose", action="store_true", help="Show extra info on terminal")
    #parser.add_argument("--isolated_nodes", action="store_true", help="Displays nodes of degree zero on the final graph")
    parser.add_argument("--isometry", action="store_true", help="Applies a random isometry to the original tensor and displays it")
    
    return parser.parse_args()

def gen_graph(T, n,m,k, F, deg_0, l_bound, u_bound,verbose):
    start = time.time()
    G = tensor_to_graph(T, n, m, k, F)
    print(f"Computation time: {time.time() - start}")
    print("Tensor T:")
    for i in range(n):
        print(T[i])
        print()
    
    #Serialize and save graph into graph.txt
    # TODO

    if verbose:
        print("\nGraph vertices: ")
        print(G.vertices())
        print("\nGraph edges: ")
        print(G.edges())

    #Identify vertices inside upper and lower bound
    print("Removing vertices of out-of-range degree")
    if deg_0:
        l_bound = -1
    out_of_bounds = []
    for v in G.vertices():
        deg = G.degree(v)
        if deg <= l_bound or deg >= u_bound:
            out_of_bounds.append(v)
    G.delete_vertices(out_of_bounds)

    return G

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
    deg_0 = False # args.isolated_nodes 

    #degree upper and lower bound filter
    u_bound = args.deg_ubound
    l_bound = args.deg_lbound

    #test if isometry flag is on
    iso = args.isometry

    #Extract cycle size and cycle ranges
    cycle_size = args.c
    strict = args.strict

    #check passed parameters
    print(n,m,k,q,labeled, verbose, u_bound, l_bound)
    
    #Create a random 3-tensor T of dimensions n x m x k over F.
    T = [[[F.random_element() for _ in range(k)] for _ in range(m)] for _ in range(n)]
    G = gen_graph(T, n,m,k, F, deg_0, l_bound, u_bound,verbose)

    #Display graph
    graph_display(G,n,m,k,q,labeled=labeled, cycle=cycle_size, strict=strict)

    
    if iso:
        #Generate element from orbit of T
        print("Applying random isometry to tensor")
        
        A = random_matrix(F, n, n, algorithm='unimodular')
        B = random_matrix(F, m, m, algorithm='unimodular')
        C = random_matrix(F, k, k, algorithm='unimodular')


        #show isometry
        if verbose:
            print("A")
            print(A)
            print("B")
            print(B)
            print("C")
            print(C)
        
        #Apply isometry: T2 = T(A,B,C)
        T2 = apply_isometry(T, A, B, C)
        
        #Generate graph and filter nodes based on cmd line arguments
        G2 = gen_graph(T2, n,m,k, F, deg_0, l_bound, u_bound,verbose)

        #Display graph
        graph_display(G2,n,m,k,q, labeled=labeled)
    