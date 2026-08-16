# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

A from-scratch hp-finite element solver in C++ (~2,500 lines), written for an undergraduate maths dissertation (Nottingham G14DIS). Solves `-Δu = f` in 1D and on 2D triangles with a hierarchical modal basis, plus a semilinear problem `-Δu + λu^(2q+1) = f` via damped Picard iteration. Dirichlet BCs only.

`dissertation_final_report.pdf` (repo root) is the authoritative spec. §5 walks through every class and the maths it implements; §5.2.2 defines the basis functions exactly. Read the relevant section before changing numerics rather than re-deriving them.

This repo is being reworked into a portfolio piece for a GitHub profile linked on a CV. Presentation quality — a working cross-platform build, a README, a clean tree — matters alongside correctness.

## Build and run

```
cmake -B build -S . && cmake --build build
cd main && ./FEM
```

Requires system Eigen (`pacman -S eigen`; 5.0.1 is what this builds against). C++17.

**The binary must be run from `main/`.** `FE_Mesh2D::constructMesh` opens mesh files by relative path, and `solution.csv` / `hp_error.csv` are written to the CWD — so CMake deliberately puts the executable in `main/`.

`.vscode/tasks.json` still holds the original Windows MSYS2 build task. It is dead on Linux; CMake replaces it.

The build is warning-clean under `-Wall -Wextra`. Keep it that way. The stricter set the old `.vscode` config used still reports ~590, almost all `-Wsign-conversion` noise.

## Testing

None. No test framework, no test files, no CI. Correctness is checked by hand: `main/main.cpp` prints a solution against a hardcoded expected vector in a comment, and L2-error convergence sweeps are written to `main/hp_error.csv`. No linter or formatter is configured; the warning flags in `.vscode/settings.json` belong to a VS Code extension, not to any build.

Benchmarks with known analytic solutions (dissertation §7) — use these to verify changes:

| Problem | `f` | `u` |
| --- | --- | --- |
| 1D Poisson | `π²sin(πx)` | `sin(πx)` |
| 1D Poisson | `1` | `x(1-x)/2` |
| 2D Poisson | `2π²sin(πx)sin(πy)` | `sin(πx)sin(πy)` |

### Known accuracy status

Measured with `h` taken from `getMeshSize()` (see below), fitted slope of `log(error)` vs `log(h)`:

| `p` | expected | 1D | 2D |
| --- | --- | --- | --- |
| 1 | 2 | 2.00 | 1.97 |
| 2 | 3 | 3.00 | 2.95 |
| 3 | 4 | 4.00 | 4.03 |

**Both dimensions now reach the theoretical rate for `p = 1, 2, 3`**, which is the range the sweep covers. The dissertation (§6) reported 2D as "converges but not at the theoretical rate" and left it unresolved; that shortfall is fixed. It had three independent causes, all since corrected: the fabricated `h` (see Configuration and inputs), an off-by-two in the edge-mode kernel index, and the missing shared-edge sign correction.

**`p ≥ 4` is still broken** — fitted rates of roughly 0.3 (`p=4`) and 1.1 (`p=5`) against expected 5 and 6. This is a separate, pre-existing defect: it was equally broken before the above fixes (1.03 and 1.26 at the time). It is *not* a basis-structure problem — bubbles vanish on all three edges and edge functions vanish on the other two edges at every degree up to 5, all to machine precision. Raising the quadrature to `nq = p+4` cuts the `p=4` error by ~2400× but leaves the rate wrong and degrades `p=3`, and the finest-mesh errors for `p = 3, 4, 5` then cluster around 5e-06, which looks like an error floor rather than a quadrature-order problem. Unresolved; start from the floor, not from the basis.

## Architecture

`FE_Solution` (façade) → `FE_Mesh` (abstract; `FE_Mesh1D` / `FE_Mesh2D`) → `Element` (abstract; `Element1D` / `Element2D`) → `PolynomialSpace` (hierarchical Lobatto/kernel basis). Assembly fills a `CSRMatrix`, which `Solver` converts to `Eigen::SparseMatrix` and solves with `SimplicialLLT`.

`FE_Solution::solve` runs a fixed pipeline: `constructMesh → allocateStiffness → assembleStiffnessMatrix → assembleLoadVector → applyBoundaryConditions → Solver::setupEigen → Solver::solveEigen`.

Numerical conventions (dissertation §5.2) — easy to break silently:

- Reference elements are `[-1, 1]` in 1D and the triangle `{-1 < ξ, η; ξ + η < 0}` in 2D. **Not** the unit triangle; affine coordinates and quadrature both assume this.
- The basis is modal/hierarchical, not nodal: vertex functions, then `p-1` Lobatto edge modes per edge (`p ≥ 2`), then `(p-1)(p-2)/2` interior bubbles (`p ≥ 3`), indexed in that order.
- Edge modes run `k = 2, ..., p`, and `evaluate_edge`/`kernel` are indexed `k-2`. The loop counter in `basis_2D`/`basis_2D_grad` is 0-based, so it must have 2 added before use. Getting this wrong is silent: `kernel(-2)` and `kernel(-1)` return `2/(1±x)`, which are not polynomials and blow up at the edge endpoints, but the guard in `evaluate_edge` zeroes them exactly at the vertices so the functions still look plausible. The invariant to test against: along edge 0 (`η = -1`), edge mode `m` must equal the 1D Lobatto function `l_{m+2}(ξ)`.
- Quadrature uses `nq = p + 1` points, exact to degree `2·nq - 1`.
- 2D load scales by `detA/4`, but 2D stiffness scales by `detA` — even though the comment directly above it (`source/Element2D.cpp:66`) states `detA/4`. Comment and code disagree, so the comment is what needs fixing: `p=1` converges at exactly the optimal rate, which it could not if the scaling were wrong.

Sharp edges:

- **Shared-edge sign correction.** `Element2D::basisSign` implements the `-1` multiplier from dissertation §5.2.1/§5.2.2. Edge modes are parameterised from the first local vertex of the edge to the second, so neighbouring triangles must agree on that direction; the canonical direction is low global vertex index to high, and when the local ordering disagrees, modes with an **odd** kernel index flip sign (even modes are symmetric under the reversal and must not be touched). `FE_Mesh2D::makeEdgePair` already sorts the pair so both triangles share DoF indices — the sign is the other half of that. It must be applied at *every* basis evaluation, not just assembly: local stiffness and load, the semilinear stiffness-product and nonlinear load, and the three reconstruction sites in `FE_Mesh2D` (`evaluateSolution`, `evaluateDerivative`, `sendSolutionToFile`). Miss one and the solve and the reported error disagree.
- `CSRMatrix::operator()(i,j)` throws `std::out_of_range` unless the slot was pre-allocated by `FE_Mesh::allocateStiffness()`. Any change to DoF connectivity must be reflected in the allocation pass first.
- `CSRMatrix::row_start` holds `noRows + 1` entries. Any loop that reads `row_start[i+1]` must stop at `size() - 1`, not `size()`. Getting this wrong reads one past the end — it trips a libstdc++ assertion at `-O0` but silently corrupts the heap in the Release build, which is how it hid for so long.
- `SimplicialLLT` is a Cholesky solver and assumes SPD. Boundary conditions are applied by zeroing rows/columns and setting a unit diagonal specifically to preserve this; non-Dirichlet BCs would break the solver choice.
- Constructor parameters shadow member names throughout (`Element`, `FE_Mesh`, `FE_Solution`). Initializer lists now match declaration order — keep them that way, since with the shadowing in place a mismatched list is easy to misread as reading a member when it reads a parameter.
- The semilinear path is 2D-only (`dynamic_cast<FE_Mesh2D*>`, unchecked) and has an unresolved `// TODO` for applying BCs inside the Picard loop.

## Configuration and inputs

Problem setup is compile-time only: `polynomialDegree`, `dimension`, and `semilinear` are literals at the top of `main/main.cpp`. There are no CLI arguments and no environment variables anywhere in the codebase.

2D meshes come from the external Triangle generator, which is not in the repo. They are selected by bare number — `"6"` loads `main/domain.6.node` and `.ele` — where the number is a refinement level (`h = 1/2^(k+1)`). Triangle files are 1-indexed and decremented on read; only 3-node triangles are supported. The Triangle invocation that produced them was not recorded. `domain.L.*` is the L-shaped domain from dissertation §7.2 — results for it exist, so it was run by editing `main.cpp`, even though no current code path names it.

**Never infer `h` from the mesh filename.** The original sweep used `h = 1/2^(k+1)`, which was wrong by 11–64× and, worse, assumed successive meshes differ by a factor of 2 when they actually differ by √2 — halving every slope on the log-log plot. Use `FE_Solution::getMeshSize()`, which measures the largest element diameter. Actual sizes: `domain.4`–`domain.8` are 32/64/128/256/512 elements at `h` = 0.354, 0.250, 0.177, 0.125, 0.088. The family alternates between two topologies, so compare meshes two apart for a clean factor-2 step. `domain.9` has 15565 elements at `h = 0.0625` — a graded mesh, not a uniform refinement of the others — and is excluded from the sweep.

## Style

Hard tabs, Allman braces, `#ifndef CLASSNAMEHEADERDEF` guards (no `#pragma once`), everything `public:`, explicit destructors with `virtual` on bases. Comments are lowercase `//` lines above blocks, usually stating the formula being implemented; no Doxygen.

Naming is inconsistent: methods are `camelCase` in the mesh/solver layer but `snake_case` in the polynomial layer, and `source/CSR_Matrix.cpp` mixes tabs with 4-space indentation. Match whatever the file you are editing already does — do not normalize across the repo for now.

## Repo etiquette

Commit messages are lowercase, free-form, descriptive sentences — no conventional-commit prefixes, scopes, or issue refs. Feature branches are `PascalCase` and merged into `main` via GitHub PRs.

## Portfolio cleanup candidates

Worth flagging when relevant: `FEM.out` is a committed 4.6 MB Windows PE32+ binary (stale, not runnable here, still in git history); there is no README; and `.vscode/` is tracked despite being gitignored.

Extensions the dissertation itself proposes (§8), in its own order of preference: Neumann and mixed boundary conditions; hp-adaptivity driven by a posteriori error estimates; fully nonlinear problems.
