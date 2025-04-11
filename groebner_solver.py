import argparse
from sage.all import *
import random
#---------------------------
# Edge condition (polynomial) helper functions.
#---------------------------
def edge_UV(C, u, v, n, m, k, R):
    """
    For l < k, computes the polynomials P_l(u,v) = sum_{i=0}^{n-1} sum_{j=0}^{m-1} C[i][j][l] * u[i] * v[j]
    u,v is an edge iff for all l < k, P_l(u,v) = 0
    """
    eqs = []
    for l in range(k):
        poly = R(0)
        for i in range(n):
            for j in range(m):
                poly += C[i][j][l] * u[i] * v[j]
        eqs.append(poly)
    return eqs

def edge_UW(C, u, w, n, m, k, R):
    """
    For each j in 0,...,m-1 computes P_j(u,w) = sum_{i=0}^{n-1} sum_{l=0}^{k-1} C[i][j][l] * u[i] * w[l]
    u,w is an edge if P_j(u,w) = 0 for all j
    """
    eqs = []
    for j in range(m):
        poly = R(0)
        for i in range(n):
            for l in range(k):
                poly += C[i][j][l] * u[i] * w[l]
        eqs.append(poly)
    return eqs

def edge_VW(C, v, w, n, m, k, R):
    """
    For each i in 0,...,n-1 computes P_i(v,w) = sum_{j=0}^{m-1} sum_{l=0}^{k-1} C[i][j][l] * v[j] * w[l]
    v,w is an edge if for all i P_i(v,w) = 0
    """
    eqs = []
    for i in range(n):
        poly = R(0)
        for j in range(m):
            for l in range(k):
                poly += C[i][j][l] * v[j] * w[l]
        eqs.append(poly)
    return eqs


#---------------------------
# Cycle type functions.
#---------------------------

def typeA_closed_walks(C, n, m, k, q):
    """
    Type A: vertices: u, u' in U; v, v' in V.

    we enforce u_0 = u'_1 = 1 and u_1 = u'_0 = 0
    and v_0 = v'_1 = 1 and v_1 = v'_0 = 0 to ensure ideal of degree 0
    """
    var_names = (['x{}'.format(i) for i in range(3, n+1)] +
                 ['y{}'.format(i) for i in range(3, n+1)] +
                 ['z{}'.format(j) for j in range(3, m+1)] +
                 ['w{}'.format(j) for j in range(3, m+1)])
    R = PolynomialRing(GF(q), var_names, order='degrevlex')
    gens = R.gens()
    u      = [R(1), R(0)] + list(gens[0:(n-2)])                   # u in P(U)
    uprime = [R(0), R(1)] + list(gens[(n-2):(2*(n-2))])            # u' in P(U)
    v      = [R(1),R(0)] + list(gens[(2*(n-2)):(2*(n-2)+(m-2))])   # v in P(V)
    vprime = [R(0), R(1)] + list(gens[(2*(n-2)+(m-2)):(2*(n-2)+2*(m-2))])  # v' in P(V)
    
    eqs = []
    eqs += edge_UV(C, u, v, n, m, k, R)            # edge u -> v
    eqs += edge_UV(C, uprime, v, n, m, k, R)         # edge v -> u'
    eqs += edge_UV(C, uprime, vprime, n, m, k, R)     # edge u' -> v'
    eqs += edge_UV(C, u, vprime, n, m, k, R)         # edge v' -> u
    
    I = R.ideal(eqs)
    I.groebner_basis()
    if verbose:
        print("Computing variety")
    sols = I.variety()
    return sols

# Type B: vertices: u, u' in U; w, w' in W.
def typeB_closed_walks(C, n, m, k, q):
    """
    Type B: vertices: u, u' in U; w, w' in W

    we enforce u_0 = u'_1 = 1 and u_1 = u'_0 = 0
    and w_0 = w'_1 = 1 and w_1 = w'_0 = 0 to ensure ideal of degree 0
    """
    var_names = (['x{}'.format(i) for i in range(3, n+1)] +
                 ['y{}'.format(i) for i in range(3, n+1)] +
                 ['z{}'.format(l) for l in range(3, k+1)] +
                 ['w{}'.format(l) for l in range(3, k+1)])
    R = PolynomialRing(GF(q), var_names, order='degrevlex')
    gens = R.gens()
    u      = [R(1), R(0)] + list(gens[0:(n-2)])                   # u in P(U)
    uprime = [R(0), R(1)] + list(gens[(n-2):(2*(n-2))])            # u' in P(U)
    w      = [R(1), R(0)] + list(gens[(2*(n-2)):(2*(n-2)+(k-2))])   # w in P(W)
    wprime = [R(0), R(1)] + list(gens[(2*(n-2)+(k-2)):(2*(n-2)+2*(k-2))])  # w' in P(W)
    
    eqs = []
    eqs += edge_UW(C, u, w, n, m, k, R)            # edge u -> w
    eqs += edge_UW(C, uprime, w, n, m, k, R)         # edge w -> u'
    eqs += edge_UW(C, uprime, wprime, n, m, k, R)     # edge u' -> w'
    eqs += edge_UW(C, u, wprime, n, m, k, R)         # edge w' -> u

    
    I = R.ideal(eqs)
    if verbose:
        print("Solving variety")
    sols = I.variety()
    return sols

# Type C: vertices: v, v' in V; w, w' in W.
def typeC_closed_walks(C, n, m, k, q):
    """
    Type C: vertices: v, v' in V; w, w' in W.

    we enforce v_0 = v'_1 = 1 and v_1 = v'_0 = 0
    and w_0 = w'_1 = 1 and w_1 = w'_0 = 0 to ensure ideal of degree 0
    """
    var_names = (['x{}'.format(j) for j in range(3, m+1)] +
                 ['y{}'.format(j) for j in range(3, m+1)] +
                 ['z{}'.format(l) for l in range(3, k+1)] +
                 ['w{}'.format(l) for l in range(3, k+1)])
    R = PolynomialRing(GF(q), var_names, order='degrevlex')
    gens = R.gens()
    v      = [R(1), R(0)] + list(gens[0:(m-2)])                   # v in P(V)
    vprime = [R(0), R(1)] + list(gens[(m-2):(2*(m-2))])            # v' in P(V)
    w      = [R(1), R(0)] + list(gens[(2*(m-2)):(2*(m-2)+(k-2))])   # w in P(W)
    wprime = [R(0), R(1)] + list(gens[(2*(m-2)+(k-2)):(2*(m-2)+2*(k-2))])  # w' in P(W)
    
    eqs = []
    eqs += edge_VW(C, v, w, n, m, k, R)            # edge v -> w
    eqs += edge_VW(C, vprime, w, n, m, k, R)         # edge w -> v'
    eqs += edge_VW(C, vprime, wprime, n, m, k, R)     # edge v' -> w'
    eqs += edge_VW(C, v, wprime, n, m, k, R)         # edge w' -> v

    I = R.ideal(eqs)
    if verbose:
        print("Solving variety")
    sols = I.variety()
    return sols

# Type D: vertices: u, u' in U; v in V; w in W.
def typeD_closed_walks(C, n, m, k, q):
    """
    Type D: vertices: u, u' in U; v in V; w in W

    we enforce the first coordinate of v and w to be 1
    and u_0 != u'_0 
    """
    var_names = (['x{}'.format(i) for i in range(3, n+1)] +
                 ['y{}'.format(i) for i in range(3, n+1)] +
                 ['z{}'.format(j) for j in range(2, m+1)] +
                 ['w{}'.format(l) for l in range(2, k+1)])
    R = PolynomialRing(GF(q), var_names, order='degrevlex')
    gens = R.gens()
    u      = [R(1), R(0)] + list(gens[0:(n-2)])
    uprime = [R(0), R(1)] + list(gens[(n-2):(2*(n-2))])
    v      = [R(1)] + list(gens[(2*(n-2)):(2*(n-2)+(m-1))])
    w      = [R(1)] + list(gens[(2*(n-2)+(m-1)):(2*(n-2)+(m-1)+(k-1))])
    
    eqs = []
    eqs += edge_UV(C, u, v, n, m, k, R)            # edge u -> v
    eqs += edge_UV(C, uprime, v, n, m, k, R)         # edge v -> u'
    eqs += edge_UW(C, uprime, w, n, m, k, R)         # edge u' -> w
    eqs += edge_UW(C, u, w, n, m, k, R)              # edge w -> u
    
    I = R.ideal(eqs)
    if verbose:
        print("Solving variety")
    sols = I.variety()
    return sols

# Type E: vertices: u in U; v, v' in V; w in W.
def typeE_closed_walks(C, n, m, k, q):
    """
    Type E: vertices: u in U; v, v' in V; w in W.

    we enforce the first coordinate of u and w to be 1
    and v_0 != v'_0 
    """
    var_names = (['x{}'.format(i) for i in range(2, n+1)] +
                 ['z{}'.format(j) for j in range(3, m+1)] +
                 ['w{}'.format(j) for j in range(3, m+1)] +
                 ['y{}'.format(l) for l in range(2, k+1)])
    R = PolynomialRing(GF(q), var_names, order='degrevlex')
    gens = R.gens()
    u      = [R(1)] + list(gens[0:(n-1)])
    v      = [R(1), R(0)] + list(gens[(n-1):(n-1)+(m-2)])
    vprime = [R(0), R(1)] + list(gens[(n-1)+(m-2):(n-1)+2*(m-2)])
    w      = [R(1)] + list(gens[(n-1)+2*(m-2):(n-1)+2*(m-2)+(k-1)])
    
    eqs = []
    eqs += edge_UV(C, u, v, n, m, k, R)            # u -> v
    eqs += edge_VW(C, v, w, n, m, k, R)            # v -> w
    eqs += edge_VW(C, vprime, w, n, m, k, R)         # w -> v'
    eqs += edge_UV(C, u, vprime, n, m, k, R)         # v' -> u
    
    I = R.ideal(eqs)
    if verbose:
        print("Solving variety")
    sols = I.variety()
    return sols

# Type F: vertices: u in U; v in V; w, w' in W.
def typeF_closed_walks(C, n, m, k, q):
    """
    Type F: vertices: u in U; v in V; w, w' in W.

    we enforce the first coordinate of u and v to be 1
    and w_0 != w'_0 
    """
    var_names = (['x{}'.format(i) for i in range(2, n+1)] +
                 ['z{}'.format(j) for j in range(2, m+1)] +
                 ['y{}'.format(l) for l in range(3, k+1)] +
                 ['w{}'.format(l) for l in range(3, k+1)])
    R = PolynomialRing(GF(q), var_names, order='degrevlex')
    gens = R.gens()
    u      = [R(1)] + list(gens[0:(n-1)])
    v      = [R(1)] + list(gens[(n-1):(n-1)+(m-1)])
    w      = [R(1), R(0)] + list(gens[(n-1)+(m-1):(n-1)+(m-1)+(k-2)])
    wprime = [R(0), R(1)] + list(gens[(n-1)+(m-1)+(k-2):(n-1)+(m-1)+2*(k-2)])
    
    eqs = []
    eqs += edge_UW(C, u, w, n, m, k, R)            # u -> w
    eqs += edge_VW(C, v, w, n, m, k, R)            # w -> v
    eqs += edge_VW(C, v, wprime, n, m, k, R)         # v -> w'
    eqs += edge_UW(C, u, wprime, n, m, k, R)         # w' -> u
    
    I = R.ideal(eqs)
    if verbose:
        print("Solving variety")
    sols = I.variety()
    return sols

#---------------------------
# 3. Master function to run all types
#---------------------------
def find_all_4cycles(C, n, m, k, q):
    sols = [0] * 6
    functions = [("A",typeA_closed_walks), \
                 ("B", typeB_closed_walks), \
                 ("C", typeC_closed_walks), \
                 ("D", typeD_closed_walks), \
                 ("E", typeE_closed_walks), \
                 ("F", typeF_closed_walks)]
    
    for idx, (type_, func) in enumerate(functions):
        if verbose:
            print(f"Computing type {type_}")
        sols[idx] = func(C,n,m,k,q)
        if verbose:
            print(f"{len(sols[idx])} walks of type {type_} found")
    
    return {"Type A": sols[0], "Type B": sols[1], "Type C": sols[2],
            "Type D": sols[3], "Type E": sols[4], "Type F": sols[5]}

#---------------------------
# 4. Usage
#---------------------------
def example_all_types(q,n,m,k):
    GFq = GF(q)
    
    random.seed(0)
    C = [[[GFq.random_element() for _ in range(k)]
             for __ in range(m)]
             for ___ in range(n)]
    
    if verbose:
        print("Random 3-tensor C:")
        for i in range(n):
            for j in range(m):
                print(C[i][j])
    
    solutions = find_all_4cycles(C, n, m, k, q)
    
    if verbose:
        print("Enlisting all solutions")
    total_sol = 0
    for typ, sol in solutions.items():
        if verbose:
            print(f"{typ}: Found {len(sol)} solution(s)")
        elif csv:
            print(f"{len(sol)},",end="") 
        if sol:
            print(f"{typ}: {sol}")
            total_sol += len(sol)
    if not(csv):
        print(f"Total number of 4-cycles: {total_sol}")
    else:
        print()

def argparser():
    #Parses the values of (n,m,k,q,labeled) as described
    parser = argparse.ArgumentParser(
        description="Finds constrained 4-cycles of random tensor graph",
        epilog="Example usage: sage groebner_solver.py -q=5 -n=4 --same_dim --minimal"
    )
    parser.add_argument("-n", type=int, default=5, help="Dimension n for the first vector space")
    parser.add_argument("-m", type=int, default=5, help="Dimension m for the second vector space")
    parser.add_argument("-k", type=int, default=5, help="Dimension k for the third vector space")
    parser.add_argument("-q", type=int, default=13, help="Prime field size")
    parser.add_argument("--same_dim", action="store_true", help="Coerces each dimension to be equal to n (i.e. n = m = k)")
    parser.add_argument("--minimal", action="store_true", help="Only displays random tensor and solutions (if any)")
    parser.add_argument("--csv", action="store_true", help="Outputs results in csv format")
    return parser.parse_args()

if __name__ == "__main__":
    global verbose 
    global csv
    args = argparser()
    
    #Dimensions of U, V, W.
    n = args.n
    m = args.m
    k = args.k
    #Field size
    q = args.q
    #equal dim
    if args.same_dim:
        m = k = n
    #verbose:
    verbose = not(args.minimal)
    #csv
    csv = args.csv
    if csv:
        print(f"{n},{m},{k},{q},",end="")
        verbose = False
    # Run the example
    example_all_types(q,n,m,k)