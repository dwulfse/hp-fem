# hp-fem

A higher-order finite element solver written from scratch in C++, for second order
elliptic boundary value problems in one and two dimensions.

Solves the Poisson problem

$$-\Delta u = f \quad \text{in } \Omega, \qquad u = 0 \quad \text{on } \partial\Omega,$$

and the semilinear problem

$$-\Delta u + \lambda u^{2q+1} = f,$$

the latter by a damped Picard iteration. Approximation uses a hierarchical modal basis
built from integrated Legendre polynomials, so the polynomial degree `p` is a parameter
rather than a rewrite. No finite element library is used; Eigen supplies only the sparse
Cholesky solve.

![L2 error under h-refinement](docs/convergence.svg)

Both dimensions attain the theoretical $O(h^{p+1})$ rate in the $L^2$ norm.

| `p` | theoretical rate | 1D measured | 2D measured |
| --- | --- | --- | --- |
| 1 | 2 | 2.00 | 1.97 |
| 2 | 3 | 3.00 | 2.95 |
| 3 | 4 | 4.00 | 4.03 |
| 4 | 5 | 5.00 | 4.97 |

Rates are least-squares fits of $\log \lVert u - u_h \rVert_{L^2}$ against $\log h$, with
`h` measured as the largest element diameter. On the finest mesh at `p = 4` the error is
$8.5 \times 10^{-10}$ in 1D and $4.4 \times 10^{-8}$ in 2D.

## Building

Requires a C++17 compiler, CMake 3.16+, and Eigen. Built and tested against Eigen 5.0.1;
no version constraint is imposed, and only the stable `Sparse` and `SparseCholesky`
interfaces are used.

```sh
cmake -B build -S .
cmake --build build
```

On Arch, Eigen comes from `pacman -S eigen`; on Debian or Ubuntu, `apt install libeigen3-dev`.

## Running

The binary is written into `main/` and must be run from there, since it reads mesh files
by relative path and writes its output to the working directory.

```sh
cd main && ./FEM
```

This writes `solution.csv` (the solution sampled on a uniform grid, alongside the
analytic solution) and `hp_error.csv` (the convergence sweep).

The problem is configured at compile time, at the top of `main/main.cpp`:

```cpp
const int polynomialDegree = 1;
const int dimension        = 1;
const bool semilinear      = false;
```

To regenerate the convergence figure from a sweep:

```sh
python3 scripts/plot_convergence.py
```

The script uses only the standard library.

## Method

**Reference elements.** The interval $[-1, 1]$ in 1D, and the triangle
$\lbrace (\xi, \eta) : -1 < \xi, \eta ; \ \xi + \eta < 0 \rbrace$ in 2D. Physical elements
are reached by an affine map, with gradients transformed by $D^{-T}$.

**Basis.** Hierarchical and modal rather than nodal, following Šolín et al. Each element
carries three vertex functions, $p - 1$ Lobatto edge modes per edge for $p \ge 2$, and
$(p-1)(p-2)/2$ interior bubbles for $p \ge 3$. Edge modes on a shared edge are
parameterised consistently by ordering the global vertex indices, with a sign correction
where the local ordering disagrees, so that the space stays conforming at odd degree.

**Quadrature.** Gauss–Legendre with $n_q = p + 1$ points, exact to degree $2n_q - 1$;
nodes above the tabulated cases are found by Newton iteration on the Legendre polynomials.
The triangle rule is a collapsed tensor product.

**Assembly and solution.** Local contributions are accumulated into a compressed sparse
row matrix whose sparsity pattern is precomputed from element connectivity. Dirichlet
conditions are imposed by zeroing the corresponding rows and columns and setting a unit
diagonal, which preserves symmetric positive definiteness, so the system is solved by the
sparse Cholesky factorisation `Eigen::SimplicialLLT`.

**Semilinear problems.** Given a relaxation parameter $0 < \alpha < 1$, the iteration

$$a(U_{n+1}, v) = (1 - \alpha)\, a(U_n, v) + \alpha \left( \int_\Omega f v \, dx - \lambda\, b(U_n; v) \right)$$

linearises the problem at each step, and the resulting system is solved directly.
Convergence is measured in the energy norm. This path is two-dimensional only.

## Layout

```
include/    class declarations
source/     implementations
main/       driver, meshes, and output
scripts/    convergence plotting
docs/       figure and the data behind it
```

The design separates the mesh, the element, and the polynomial space:
`FE_Solution` drives the solve and owns an `FE_Mesh` (`FE_Mesh1D` or `FE_Mesh2D`), which
holds `Element` objects (`Element1D` or `Element2D`), each of which owns a
`PolynomialSpace`. Dimension-specific behaviour lives behind the abstract bases, so the
assembly and solve paths are shared.

## Meshes

Two-dimensional meshes are generated externally by
[Triangle](https://www.cs.cmu.edu/~quake/triangle.html) and read from its `.node` and
`.ele` files. They are selected by number in `main/main.cpp` — `"6"` loads
`main/domain.6.node` and `main/domain.6.ele`. Only three-node triangles are supported.
`main/domain.L.*` is an L-shaped domain.

## Background

Written for a fourth-year mathematics dissertation at the University of Nottingham
(G14DIS), supervised by Dr. Paul Houston. `dissertation_final_report.pdf` in this
repository gives the full derivation: function spaces and weak formulations in §3–4, the
implementation and every basis function in §5, and convergence studies in §6.

Principal references are Šolín, Segeth and Doležel, *Higher-Order Finite Element Methods*,
for the hierarchical triangular elements; Karniadakis and Sherwin, *Spectral/hp Element
Methods for Computational Fluid Dynamics*; and Ciarlet, *The Finite Element Method for
Elliptic Problems*, for the underlying theory.

## Scope

Dirichlet boundary conditions only. Meshes are read, not generated or adapted. The natural
extensions, in order of appeal, are Neumann and mixed conditions, `hp`-adaptivity driven by
a posteriori error estimates, and fully nonlinear problems.
