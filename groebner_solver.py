#=====================================================================
# SageMath script: Finding all 4–cycles in the 3–tensor graph via Groebner basis
# with full projective normalization.
#
# We work in projective coordinates by representing a point in P(U) as (1,x2,...,xn),
# in P(V) as (1,z2,...,zm), and in P(W) as (1,w2,...,wk).
#
# For each cycle type, extra linear constraints are added to remove any residual free
# parameters so that the ideal becomes 0-dimensional.
#=====================================================================
import argparse
from sage.all import *

#---------------------------
# 1. Edge condition helper functions.
#---------------------------
def edge_UV(C, u, v, n, m, k, R):
    """
    For vertices u in U (length n) and v in V (length m), returns a list of k polynomials in R enforcing:
      For each l in 0,...,k-1:
         sum_{i=0}^{n-1} sum_{j=0}^{m-1} C[i][j][l] * u[i] * v[j] = 0.
    u and v are given in projective coordinates (first coordinate = 1).
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
    For vertices u in U (length n) and w in W (length k), returns a list of m polynomials in R enforcing:
      For each j in 0,...,m-1:
         sum_{i=0}^{n-1} sum_{l=0}^{k-1} C[i][j][l] * u[i] * w[l] = 0.
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
    For vertices v in V (length m) and w in W (length k), returns a list of n polynomials in R enforcing:
      For each i in 0,...,n-1:
         sum_{j=0}^{m-1} sum_{l=0}^{k-1} C[i][j][l] * v[j] * w[l] = 0.
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
# 2. Cycle type functions.
# Each vertex is represented with its first coordinate fixed to 1.
# Extra linear constraints are added to force the system to be 0-dimensional.
#---------------------------

# Type A: vertices: u, u' in U; v, v' in V.
def typeA_closed_walks(C, n, m, k, q):
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
    #assert (I.dimension() == 0)
    #print("Computing Groebner basis using F4 for Type A")
    #I.groebner_basis(algorithm='slimgb')
    if verbose:
        print("Groebner basis computed. Now computing variety.")
    sols = I.variety()
    return sols

# Type B: vertices: u, u' in U; w, w' in W.
def typeB_closed_walks(C, n, m, k, q):
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
# 3. Master function to run all types and display solution counts.
#---------------------------
def find_all_4cycles(C, n, m, k, q):
    if verbose:
        print("Computing type A")
    solA = typeA_closed_walks(C, n, m, k, q)
    if verbose:
        print(f"{len(solA)} walks of type A found")
        print("Computing type B")
    solB = typeB_closed_walks(C, n, m, k, q)
    if verbose:
        print(f"{len(solB)} walks of type B found")
        print("Computing type C")
    solC = typeC_closed_walks(C, n, m, k, q)
    if verbose:
        print(f"{len(solC)} walks of type C found")
        print("Computing type D")
    solD = typeD_closed_walks(C, n, m, k, q)
    if verbose:
        print(f"{len(solD)} walks of type D found")
        print("Computing type E")
    solE = typeE_closed_walks(C, n, m, k, q)
    if verbose:
        print(f"{len(solE)} walks of type E found")
        print("Computing type F")
    solF = typeF_closed_walks(C, n, m, k, q)
    if verbose:
        print(f"{len(solF)} walks of type F found")
    
    return {"Type A": solA, "Type B": solB, "Type C": solC,
            "Type D": solD, "Type E": solE, "Type F": solF}

#---------------------------
# 4. Usage
#---------------------------
def example_all_types(q,n,m,k):
    GFq = GF(q)
    
    import random
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
        if sol:
            print(f"{typ}: {sol}")
            total_sol += len(sol)
    print(f"Total number of 4-cycles: {total_sol}")

def argparser():
    #Parses the values of (n,m,k,q,labeled) as described
    parser = argparse.ArgumentParser(
        description="Build graph from a random 3-tensor",
        epilog="Example usage: sage script.sage -n=4 -m=3 -k=5 -q=7 -c=3 --strict --labeled --verbose"
    )
    parser.add_argument("-n", type=int, default=5, help="Dimension n for the first vector space")
    parser.add_argument("-m", type=int, default=5, help="Dimension m for the second vector space")
    parser.add_argument("-k", type=int, default=5, help="Dimension k for the third vector space")
    parser.add_argument("-q", type=int, default=13, help="Prime field size")
    parser.add_argument("--same_dim", action="store_true", help="Coerces each dimension to be equal to n (i.e. n = m = k)")
    parser.add_argument("--minimal", action="store_true", help="Only displays random tensor and solutions (if any)")
    return parser.parse_args()

if __name__ == "__main__":
    global verbose 

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
    # Run the example
    example_all_types(q,n,m,k)