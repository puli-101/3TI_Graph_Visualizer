# 3-Tensor Graph Visualizer

## 0. Description

<b>Summary: </b>Generates a random 3-tensor represented as a 3-dimensional array, computes the associated graph and displays the result after applying a set of filters (see section 2.) 

In particular, given the finite dimensional vector spaces $U,V,W$ of dimensions $n,m,k$ over $\mathbb{F}_q$, we represent a 3-tensor as the trilinear form $\mathcal{C} : U \times V \times W \rightarrow \mathbb{F}_q$. 

The computed set of vertices $\mathcal{V}(\mathcal{C})$, corresponds to the disjoint union of projective spaces $\mathcal{V}(\mathcal{C}) := \mathbb{P}(U) \cup \mathbb{P}(V) \cup \mathbb{P}(W)$. On a high level, the set of edges $\mathcal{E}(\mathcal{C})$ is computed as the set of pairs of representatives $(u,v) \in \mathbb{P}(U) \times \mathbb{P}(V)$ (for instance), for which the missing coordinate vanishes, e.g. $(u,v) \in \mathcal{E}(\mathcal{C}) \iff \forall w \in W \mathcal{C}(u,v,w) = 0$ (for implementation details see section 4.)

## 1. Requirements

- [SageMath](https://www.sagemath.org/)

## 2. Command-line arguments


| Option | Description | Default |
|:-------|:------------|:--------|
| -n  N  | Dimension n for the first vector space |  4|
| -m  M  | Dimension m for the second vector space |  4| 
| -k  K  | Dimension k for the third vector space | 4 |
| -q  Q  | Prime field size | 5 |
|-c C    | Highlights all cycles of length c in the final graph | None |
|--loose | Highlights all cycles of length 2 < c' <= c | false |
| --deg_lbound D | Filters out all nodes of degree less or equal than specified | 0 |
| --deg_ubound D | Filters out all nodes of degree greater or equal than specified | 1000 |
| --isolated_nodes | Displays nodes of degree zero on the final graph | false |
| --labeled  | Show graph with vertex labels | false |
| --verbose | Show extra info on terminal | false | 
|--isometry | Applies a random isometry to the original tensor and displays it | false |
|--no_visualization | Prevents the display of a new tab with the resulting tensor graph | False |
|--minimal | Only prints in terminal the random tensor, and, if enabled, all cycles of various lengths | false | 
|--load_tensor | Loads tensor from file instead of generating a random one [TODO]| "" |
|--load_graph | Loads tensor graph from file instead of calculating one given a random tensor [TODO]| "" |

## 3. Sample execution

    sage main.py -n=5 -m=5 -k=5 -q=7 --deg_lbound=1 --deg_ubound=10 --verbose

## 4. Implementation and Specifications 

### 4.1 Vanishing test

Given $(u,v) \in \mathbb{P}(U) \times \mathbb{P}(V)$, the vanishing test for this couple of points is the following: 

If $\forall l \in [1,k], \sum_{i = 1}^{n} \sum_{j = 1}^{m} u_i \cdot v_j \cdot \mathcal{C}_{i,j,l} = 0$, then $(u,v) \in \mathcal{E}(\mathcal{C})$


### 4.2 Limitations

This implementation is only suited for small values $n,m,k,q$. No tests have been done with values higher than 10.

## Authors

Developed by [David Pulido Cornejo](https://github.com/puli-101) under the supervision of [Laurane Marco](https://lauranemarco.github.io/) as part of a combinatorial-algebraic study of the 3-Tensor Isomorphism Problem @ [EPFL/LASEC](https://lasec.epfl.ch/).