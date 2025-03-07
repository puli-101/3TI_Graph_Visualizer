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
def tensor_to_graph(T, n, m, k, F):
    #List elements of the projective spaces for U, V, and W.
    P_U = list(ProjectiveSpace(n-1, F))
    P_V = list(ProjectiveSpace(m-1, F))
    P_W = list(ProjectiveSpace(k-1, F))
    
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
    
    print("Labeling all vertices")
    for idx, p in enumerate(P_U):
        label = ("U", idx)
        vertices_u.append(label)
        vv_map[label] = vector(F, p)  
    for idx, p in enumerate(P_V):
        label = ("V", idx)
        vertices_v.append(label)
        vv_map[label] = vector(F, p)
    for idx, p in enumerate(P_W):
        label = ("W", idx)
        vertices_w.append(label)
        vv_map[label] = vector(F, p)
    
    vertices = vertices_u + vertices_v + vertices_w
    #Create an graph and add all vertices
    G = Graph(multiedges=False)
    G.add_vertices(vertices)
    
    #Add edges
    print("Adding U V edges")
    #Edge between a vertex from U and one from V if C(u,v,-) = 0
    for u_label in vertices_u:
        for v_label in vertices_v:
            if is_edge_UV(T, vv_map[u_label], vv_map[v_label]):
                G.add_edge(u_label, v_label)
    print("Adding U W edges")
    #Edge between U and W 
    for u_label in vertices_u:
        for w_label in vertices_w:
            if is_edge_UW(T, vv_map[u_label], vv_map[w_label]):
                G.add_edge(u_label, w_label)
    print("Adding V W edges")
    #Edge between V and W 
    for v_label in vertices_v:
        for w_label in vertices_w:
            if is_edge_VW(T, vv_map[v_label], vv_map[w_label]):
                G.add_edge(v_label, w_label)
    
    return G