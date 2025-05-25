# Closed Walk Finder in Tensor Graphs

This script identifies all **normalized closed walks of length 4** in a tripartite tensor graph derived from a random 3-tensor over a finite field 𝔽_q. It does so by algebraically encoding edge conditions between projective points and solving the resulting polynomial systems using Gröbner bases.

## Requirements

- [SageMath](https://www.sagemath.org/)

## Usage

```bash
sage groebner_solver.py [options]
```

### Example

```bash
sage groebner_solver.py -q=7 -n=4 --same_dim
```

### Arguments

| Option         | Description                                                  |
|----------------|--------------------------------------------------------------|
| `-n`           | Dimension of the first space **U** (default: 5)              |
| `-m`           | Dimension of the second space **V** (default: 5)             |
| `-k`           | Dimension of the third space **W** (default: 5)              |
| `-q`           | Prime field size **q** (default: 13)                         |
| `--same_dim`   | Sets `n = m = k`                                             |
| `--minimal`    | Only displays the random tensor and the number of solutions |
| `--csv`        | Outputs solution counts in a CSV row format                  |

## Walk Types

The solver searches for 4-cycles of the following types:

- **Type A**: `U -> V -> U′ -> V′ -> U`
- **Type B**: `U -> W -> U′ -> W′ -> U`
- **Type C**: `V -> W -> V′ -> W′ -> V`
- **Type D**: `U -> V -> U′ -> W -> U`
- **Type E**: `U -> V -> W -> V′ -> U`
- **Type F**: `U -> W -> V -> W′ -> U`

Each type encodes a different shape of closed walk with two or more pinned vertices to ensure that the ideal defined by the polynomial system is zero-dimensional (i.e., has finitely many solutions).

## Output

- By default: prints each type's walk count and solutions (if any).
- With `--minimal`: prints only the tensor and total number of solutions.
- With `--csv`: prints a CSV line like `n,m,k,q,typeA_count,typeB_count,...`.

## Warning

To use the `groebner_tester.sh` script on your system, modify the location of your sage environment

## Authors

Developed by David Pulido Cornejo under the supervision of Laurane Marco as part of a combinatorial-algebraic study of the 3-Tensor Isomorphism Problem @ EPFL/LASEC.

