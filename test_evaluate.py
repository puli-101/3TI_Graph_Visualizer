from tensor import *
from tools import *

#recursively get all permutations of a list
def get_permutations(lst):
    #base case
    if len(lst) <= 1:
        return [lst]
    
    perms = []
    #loop through each element and get permutations of the rest of the list.
    for i in range(len(lst)):
        # Current element
        current = lst[i]
        #remaining list without the current element
        remaining = lst[:i] + lst[i+1:]
        for p in get_permutations(remaining):
            perms.append([current] + p)
    return perms

if __name__ == "__main__":
    q = 5
    F = GF(q)

    C1 = [[[2, 4, 0], [1, 3, 1], [4, 4, 1]], \
        [[1, 4, 4], [1, 3, 0], [3, 0, 1]], \
        [[0, 4, 4], [1, 3, 1], [3, 3, 1]]]

    C2 = [[[2, 4, 0], [4, 2, 3], [3, 1, 0]], \
        [[3, 0, 0], [4, 3, 1], [4, 3, 1]], \
        [[1, 3, 3], [4, 2, 4], [4, 2, 2]]]

    basis = [[1,0,0],[0,1,0],[0,0,1]]

    T1 = [[4, 4, 1], [3, 1, 0], [2, 1, 1]]
    T1 = [vector(coerce_list(v,q)) for v in T1]
    T2 = [[2, 0, 1], [1, 0, 1], [4, 2, 1]]
    T2 = [vector(coerce_list(v,q)) for v in T2]

    test_triangle(coerce_tensor(C1,q),T1,basis)
    test_triangle(coerce_tensor(C2,q),T2,basis)
    A = Matrix(F, [[1,3,0], [2,1,3], [0,3,3]])
    B = Matrix(F, [[3,1,4], [4,0,0], [4,0,1]])
    C = Matrix(F, [[0,1,0], [4,0,0], [2,1,1]])

    print(f"u = {T1[0]} / A * u' = {A * T2[0]}") 
    print(f"v = {T1[1]} / B * w' = {B * T2[1]}") 
    print(f"w = {T1[2]} / C * u' = {C * T2[2]}")

