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
| 1 | 2 | 2.00 | **1.97** |
| 2 | 3 | 3.00 | 2.08 |
| 3 | 4 | 4.00 | 1.72 |

1D is optimal at every degree. **In 2D only `p=1` is optimal; `p=2` and `p=3` both stall at ~2**, i.e. the edge and bubble modes contribute nothing asymptotically. That is the signature of a non-conforming higher-order basis, and it points hard at the missing shared-edge sign correction under Sharp edges.

Note `p=1` hitting 1.97 means the vertex basis, the affine map, quadrature, assembly, BCs and the solver are all sound — so the defect is specifically in the `p ≥ 2` modes, not anywhere in the surrounding machinery.

The dissertation (§6) reports this as "converges but not at the theoretical rate". Its plot was additionally distorted by the `h` bug described under Configuration and inputs; the numbers above are post-fix and are the baseline to beat.

## Architecture

`FE_Solution` (façade) → `FE_Mesh` (abstract; `FE_Mesh1D` / `FE_Mesh2D`) → `Element` (abstract; `Element1D` / `Element2D`) → `PolynomialSpace` (hierarchical Lobatto/kernel basis). Assembly fills a `CSRMatrix`, which `Solver` converts to `Eigen::SparseMatrix` and solves with `SimplicialLLT`.

`FE_Solution::solve` runs a fixed pipeline: `constructMesh → allocateStiffness → assembleStiffnessMatrix → assembleLoadVector → applyBoundaryConditions → Solver::setupEigen → Solver::solveEigen`.

Numerical conventions (dissertation §5.2) — easy to break silently:

- Reference elements are `[-1, 1]` in 1D and the triangle `{-1 < ξ, η; ξ + η < 0}` in 2D. **Not** the unit triangle; affine coordinates and quadrature both assume this.
- The basis is modal/hierarchical, not nodal: vertex functions, then `p-1` Lobatto edge modes per edge (`p ≥ 2`), then `(p-1)(p-2)/2` interior bubbles (`p ≥ 3`), indexed in that order.
- Quadrature uses `nq = p + 1` points, exact to degree `2·nq - 1`.
- 2D load scales by `detA/4`, but 2D stiffness scales by `detA` — even though the comment directly above it (`source/Element2D.cpp:66`) states `detA/4`. Comment and code disagree, so the comment is what needs fixing: `p=1` converges at exactly the optimal rate, which it could not if the scaling were wrong.

Sharp edges:

- **Missing shared-edge sign correction — prime suspect for the 2D shortfall.** Dissertation §5.2.1/§5.2.2 prescribes a `-1` multiplier on edge basis functions when a shared edge's global vertex IDs are not in increasing order, without which odd-degree edge modes don't match across neighbouring triangles. `FE_Mesh2D::makeEdgePair` sorts the pair so both triangles get the same DoF *indices*, but no sign correction exists anywhere in the source — so the basis is likely non-conforming at odd `p`. Verify before relying on this.
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
