from sage.all import *
from sage.graphs.graph import Graph

"""
Defines tensor operations and constructs graph
"""

#Evaluate the tensor on three vectors.
def tensor_value(T, u, v, w):
    #T is assumed to be a 3-tensor given as a 3D list: T[i][j][k]
    #u, v, w are vectors 
    s = T[0][0][0].parent()(0) #zero of the field where T is defined
    for i in range(len(u)):
        for j in range(len(v)):
            for k in range(len(w)):
                s += u[i] * v[j] * w[k] * T[i][j][k]
    return s

#Test whether the partial evaluation vanishes.
def is_edge_UV(T, u, v):
    #For fixed u in U and v in V, check that
    #for every coordinate k, sum_{i,j} u_i * v_j * T[i][j][k] == 0.
    kdim = len(T[0][0])
    F = T[0][0][0].parent()
    for k in range(kdim):
        s = F(0)
        for i in range(len(u)):
            for j in range(len(v)):
                s += u[i] * v[j] * T[i][j][k]
        if s != 0:
            return False
    return True

def is_edge_UW(T, u, w):
    #For fixed u in U and w in W, check that for each j,
    #sum_{i,k} u_i * w_k * T[i][j][k] == 0.
    mdim = len(T[0])
    F = T[0][0][0].parent()
    for j in range(mdim):
        s = F(0)
        for i in range(len(u)):
            for k in range(len(w)):
                s += u[i] * w[k] * T[i][j][k]
        if s != 0:
            return False
    return True

def is_edge_VW(T, v, w):
    #For fixed v in V and w in W, check that for each i,
    #sum_{j,k} v[j] * w[k] * T[i][j][k] == 0.
    ndim = len(T)
    F = T[0][0][0].parent()
    for i in range(ndim):
        s = F(0)
        for j in range(len(v)):
            for k in range(len(w)):
                s += v[j] * w[k] * T[i][j][k]
        if s != 0:
            return False
    return True

# Main function that builds the graph associated with a 3-tensor.
def tensor_to_graph(T, n, m, k, F, verbose=False, minimal=False):
    #List elements of the projective spaces for U, V, and W.
    P_U = list(ProjectiveSpace(n-1, F))
    P_V = list(ProjectiveSpace(m-1, F))
    P_W = list(ProjectiveSpace(k-1, F))
    
    if not(minimal):
        print("Sizes of projective Spaces:")
        print(f"U : {len(P_U)}")
        print(f"V : {len(P_V)}")
        print(f"W : {len(P_W)}")

    #Create vertices for each projective point
    vertices = []
    vertices_u = []
    vertices_v = []
    vertices_w = []
    #Maps vertex label to its vector
    vv_map = {}  
    
    if not(minimal):
        print("Labeling all vertices")
    #progressive labeling with unique integer ID
    #previous versions set label as ("U",idx)/("V",idx)/("W",idx)
    id_c = 0
    for idx, p in enumerate(P_U):
        id_c += 1
        label = id_c 
        vertices_u.append(id_c)
        vv_map[label] = vector(F, p)  
    for idx, p in enumerate(P_V):
        id_c += 1
        label = id_c 
        vertices_v.append(label)
        vv_map[label] = vector(F, p)
    for idx, p in enumerate(P_W):
        id_c += 1
        label = id_c
        vertices_w.append(label)
        vv_map[label] = vector(F, p)
    
    #display label-node mapping
    if verbose:
        print("Label-node mapping")
        for k in vv_map.keys():
            print(k, vv_map[k])

    vertices = vertices_u + vertices_v + vertices_w
    #Create an graph and add all vertices
    G = Graph(multiedges=False)
    G.add_vertices(vertices)
    
    #Add edges
    if not(minimal):
        print("Adding U V edges")
    #Edge between a vertex from U and one from V if C(u,v,-) = 0
    for u_label in vertices_u:
        for v_label in vertices_v:
            if is_edge_UV(T, vv_map[u_label], vv_map[v_label]):
                G.add_edge(u_label, v_label)
    if not(minimal):
        print("Adding U W edges")
    #Edge between U and W 
    for u_label in vertices_u:
        for w_label in vertices_w:
            if is_edge_UW(T, vv_map[u_label], vv_map[w_label]):
                G.add_edge(u_label, w_label)
    if not(minimal):
        print("Adding V W edges")
    #Edge between V and W 
    for v_label in vertices_v:
        for w_label in vertices_w:
            if is_edge_VW(T, vv_map[v_label], vv_map[w_label]):
                G.add_edge(v_label, w_label)
    
    return G

def apply_isometry(T, A, B, C):
    #T : 3-tensor represented as a 3d list
    #A,B,C Invertible matrices
    #returns: T' s.t. T'(u,v,w) = T(Au,Bv,Cw) 
    #w coefficients T'_{p,q,r} = sum_{i,j,k} T_{i,j,k} * A[i,p] * B[j,q] * C[k,r]
    n = len(T)
    m = len(T[0])
    k = len(T[0][0])
    
    # Initialize the new tensor as a 3d list with zeros.
    T_prime = [[[0 for r in range(k)] for q in range(m)] for p in range(n)]
    
    # Loop over the new indices p, q, r.
    for p in range(n):
        for q in range(m):
            for r in range(k):
                s = 0
                # Sum over the old indices i, j, k_idx.
                for i in range(n):
                    for j in range(m):
                        for k_idx in range(k):
                            s += T[i][j][k_idx] * A[i, p] * B[j, q] * C[k_idx, r]
                T_prime[p][q][r] = s
    return T_prime


#========= ADDITIONAL FUNCTIONS FOR TESTING ================
#Prints evaluation of tensor C for a candidate triangle T
#tests for each element in basis if the evaluation equals 0
#Supposition: n = m = k
def test_triangle(C, T, basis):
    print("Evaluating Tensor Triangle")
    print(T)
    print("Evaluating (u,v)")
    for v in basis:
        print(f"{v}:\t{tensor_value(C,T[0],T[1], v)}")
    print(f"is_edge result: {is_edge_UV(C,T[0],T[1])}")
    print("Evaluating (u,w)")
    for v in basis:
        print(f"{v}:\t{tensor_value(C,T[0],v, T[2])}")
    print(f"is_edge result: {is_edge_UW(C,T[0],T[2])}")
    print("Evaluating (v,w)")
    for v in basis:
        print(f"{v}:\t{tensor_value(C,v,T[1], T[2])}")
    print(f"is_edge result: {is_edge_VW(C,T[1],T[2])}\n")

#Transforms 3d-array of int's to 3d-array of elements in GF(q) 
def coerce_tensor(C, q):
    return [[coerce_list(l,q) for l in M] for M in C]

#Transforms 3d-array of int's to 3d-array of elements in GF(q) 
def coerce_list(v, q):
    F = GF(q)
    return [F(x) for x in v]