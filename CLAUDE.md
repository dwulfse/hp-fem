# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

A from-scratch hp-finite element solver in C++ (~2,500 lines), written for an undergraduate maths dissertation (Nottingham G14DIS). Solves `-Δu = f` in 1D and on 2D triangles with a hierarchical modal basis, plus a semilinear problem `-Δu + λu^(2q+1) = f` via damped Picard iteration. Dirichlet BCs only.

`dissertation_final_report.pdf` (repo root) is the authoritative spec. §5 walks through every class and the maths it implements; §5.2.2 defines the basis functions exactly. Read the relevant section before changing numerics rather than re-deriving them.

This repo is being reworked into a portfolio piece for a GitHub profile linked on a CV. Presentation quality — a working cross-platform build, a README, a clean tree — matters alongside correctness.

## Build and run

```
cmake -B build -S . && cmake --build build
cd main && ./FEM --help
```

`source/` builds into a `fem_core` static library that both `FEM` and `fem_tests` link against; only `main/main.cpp` is outside it.

Requires system Eigen (`pacman -S eigen`; 5.0.1 is what this builds against). C++17.

**The binary must be run from `main/`.** `FE_Mesh2D::constructMesh` opens mesh files by relative path, and `solution.csv` / `hp_error.csv` are written to the CWD — so CMake deliberately puts the executable in `main/`.

`.vscode/` is untracked (gitignored) but still present on disk; its `tasks.json` holds the original Windows MSYS2 build task, which is dead on Linux. CMake replaces it.

The build is warning-clean under `-Wall -Wextra`. Keep it that way. The stricter set the old `.vscode` config used still reports ~590, almost all `-Wsign-conversion` noise.

## Testing

Catch2, fetched by CMake at configure time. `ctest --test-dir build --output-on-failure`, or run `main/fem_tests` directly for Catch2's own output and filtering. `-DFEM_BUILD_TESTS=OFF` skips the fetch. The binary lands in `main/` because the 2D cases read meshes by relative path.

`tests/` covers three layers: quadrature against exact integrals of monomials, the basis structurally (interpolation, vanishing, Lobatto reduction along an edge, gradients against finite differences), and the solver end to end asserting measured convergence rates.

**Every one of these was written because it catches a defect this repo actually shipped.** Before changing numerics, check whether an existing test already pins the property — and if you fix something, add the test that would have caught it. Reintroducing the quadrature typo fails four cases; reintroducing the edge kernel index fails three.

CI (`.github/workflows/ci.yml`) builds with gcc and clang in Debug and Release, runs ctest, runs the solver, and separately builds with `-Werror`. No linter or formatter is configured.

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

Verified out to `p = 5` in 2D as well (4.97 and 6.01 against expected 5 and 6), though the committed sweep only runs to `p = 3`.

**Both dimensions now reach the theoretical rate at every degree tested.** The dissertation (§6) reported 2D as "converges but not at the theoretical rate" and left it unresolved; that shortfall is fixed. It had four independent causes, all since corrected: the fabricated `h` (see Configuration and inputs), an off-by-two in the edge-mode kernel index, the missing shared-edge sign correction, and two bugs in `GaussQuadrature` (below).

When a convergence rate looks wrong, check the quadrature before the basis — the rules were the single largest source of error here, and a wrong rule degrades results without ever looking obviously broken.

## Architecture

`FE_Solution` (façade) → `FE_Mesh` (abstract; `FE_Mesh1D` / `FE_Mesh2D`) → `Element` (abstract; `Element1D` / `Element2D`) → `PolynomialSpace` (hierarchical Lobatto/kernel basis). Assembly fills a `CSRMatrix`, which `Solver` converts to `Eigen::SparseMatrix` and solves with `SimplicialLLT`.

`FE_Solution::solve` runs a fixed pipeline: `constructMesh → allocateStiffness → assembleStiffnessMatrix → assembleLoadVector → applyBoundaryConditions → Solver::setupEigen → Solver::solveEigen`.

Numerical conventions (dissertation §5.2) — easy to break silently:

- Reference elements are `[-1, 1]` in 1D and the triangle `{-1 < ξ, η; ξ + η < 0}` in 2D. **Not** the unit triangle; affine coordinates and quadrature both assume this.
- The basis is modal/hierarchical, not nodal: vertex functions, then `p-1` Lobatto edge modes per edge (`p ≥ 2`), then `(p-1)(p-2)/2` interior bubbles (`p ≥ 3`), indexed in that order.
- Edge modes run `k = 2, ..., p`, and `evaluate_edge`/`kernel` are indexed `k-2`. The loop counter in `basis_2D`/`basis_2D_grad` is 0-based, so it must have 2 added before use. Getting this wrong is silent: `kernel(-2)` and `kernel(-1)` return `2/(1±x)`, which are not polynomials and blow up at the edge endpoints, but the guard in `evaluate_edge` zeroes them exactly at the vertices so the functions still look plausible. The invariant to test against: along edge 0 (`η = -1`), edge mode `m` must equal the 1D Lobatto function `l_{m+2}(ξ)`.
- Quadrature uses `nq = p + 1` points, exact to degree `2·nq - 1`. `GaussQuadrature` hardcodes `n = 1..5` and finds Legendre roots by Newton beyond that; both paths had bugs that silently wrecked `p ≥ 4` (see Sharp edges). If you touch either, re-check the rules integrate monomials exactly — an `n`-point rule must reproduce `∫x^k` on `[-1,1]` for all `k ≤ 2n-1`. That test catches in seconds what shows up as a mystifying convergence plot otherwise.
- 2D load scales by `detA/4`, but 2D stiffness scales by `detA` — even though the comment directly above it (`source/Element2D.cpp:66`) states `detA/4`. Comment and code disagree, so the comment is what needs fixing: `p=1` converges at exactly the optimal rate, which it could not if the scaling were wrong.

Sharp edges:

- **Shared-edge sign correction.** `Element2D::basisSign` implements the `-1` multiplier from dissertation §5.2.1/§5.2.2. Edge modes are parameterised from the first local vertex of the edge to the second, so neighbouring triangles must agree on that direction; the canonical direction is low global vertex index to high, and when the local ordering disagrees, modes with an **odd** kernel index flip sign (even modes are symmetric under the reversal and must not be touched). `FE_Mesh2D::makeEdgePair` already sorts the pair so both triangles share DoF indices — the sign is the other half of that. It must be applied at *every* basis evaluation, not just assembly: local stiffness and load, the semilinear stiffness-product and nonlinear load, and the three reconstruction sites in `FE_Mesh2D` (`evaluateSolution`, `evaluateDerivative`, `sendSolutionToFile`). Miss one and the solve and the reported error disagree.
- `CSRMatrix::operator()(i,j)` throws `std::out_of_range` unless the slot was pre-allocated by `FE_Mesh::allocateStiffness()`. Any change to DoF connectivity must be reflected in the allocation pass first.
- **`GaussQuadrature` had two independent bugs, both fixed, both worth knowing about.** The hardcoded 5-point node table used `sqrt(15.0/7)` where the correct nodes are `sqrt(5/9 -+ (2/9)·sqrt(10.0/7))` — one wrong digit, and since the *weights* were right the rule still integrated constants exactly. And the Newton loop tested `abs(x - x_old)`, which resolves to the integer `abs` on a double, truncating every step below 1.0 to zero and stopping after a single iteration, so all `n ≥ 6` rules were accurate to ~1e-4 rather than 1e-14. Together these capped 2D at `p ≤ 3` (`nq = p+1`, so `p=4` is the first degree to hit the broken 5-point rule).
- `CSRMatrix::row_start` holds `noRows + 1` entries. Any loop that reads `row_start[i+1]` must stop at `size() - 1`, not `size()`. Getting this wrong reads one past the end — it trips a libstdc++ assertion at `-O0` but silently corrupts the heap in the Release build, which is how it hid for so long.
- `SimplicialLLT` is a Cholesky solver and assumes SPD. Boundary conditions are applied by zeroing rows/columns and setting a unit diagonal specifically to preserve this; non-Dirichlet BCs would break the solver choice.
- Constructor parameters shadow member names throughout (`Element`, `FE_Mesh`, `FE_Solution`). Initializer lists now match declaration order — keep them that way, since with the shadowing in place a mismatched list is easy to misread as reading a member when it reads a parameter.
- The semilinear path is 2D-only (`dynamic_cast<FE_Mesh2D*>`, unchecked) and has an unresolved `// TODO` for applying BCs inside the Picard loop.

## Configuration and inputs

Problem setup is entirely command line: `./FEM --help`. Dimension, degree, element count, mesh, problem, the semilinear parameters, `--sweep` and `--field` are all flags, parsed in `main/main.cpp`. There are no environment variables and no config files. `main.cpp` is a driver only — it holds the benchmark forcing functions, argument parsing and validation, and nothing numerical.

2D meshes come from the external Triangle generator, which is not in the repo. They are selected by bare number — `"6"` loads `main/domain.6.node` and `.ele` — where the number is a refinement level (`h = 1/2^(k+1)`). Triangle files are 1-indexed and decremented on read; only 3-node triangles are supported. The Triangle invocation that produced them was not recorded. `domain.L.*` is the L-shaped domain from dissertation §7.2 — results for it exist, so it was run by editing `main.cpp`, even though no current code path names it.

**Never infer `h` from the mesh filename.** The original sweep used `h = 1/2^(k+1)`, which was wrong by 11–64× and, worse, assumed successive meshes differ by a factor of 2 when they actually differ by √2 — halving every slope on the log-log plot. Use `FE_Solution::getMeshSize()`, which measures the largest element diameter. Actual sizes: `domain.4`–`domain.8` are 32/64/128/256/512 elements at `h` = 0.354, 0.250, 0.177, 0.125, 0.088. The family alternates between two topologies, so compare meshes two apart for a clean factor-2 step. `domain.9` has 15565 elements at `h = 0.0625` — a graded mesh, not a uniform refinement of the others — and is excluded from the sweep.

## Style

Hard tabs, Allman braces, `#ifndef CLASSNAMEHEADERDEF` guards (no `#pragma once`), everything `public:`, explicit destructors with `virtual` on bases. Comments are lowercase `//` lines above blocks, usually stating the formula being implemented; no Doxygen.

Naming is inconsistent: methods are `camelCase` in the mesh/solver layer but `snake_case` in the polynomial layer, and `source/CSR_Matrix.cpp` mixes tabs with 4-space indentation. Match whatever the file you are editing already does — do not normalize across the repo for now.

## Repo etiquette

Commit messages are lowercase, free-form, descriptive sentences — no conventional-commit prefixes, scopes, or issue refs. Feature branches are `PascalCase` and merged into `main` via GitHub PRs.

## Portfolio cleanup candidates

Mostly done: the repo is `dwulfse/hp-fem`, with a README, both figures in `docs/`, a CMake build, tests, CI, an MIT licence, and `FEM.out` / `.vscode/` / the generated CSVs all untracked.

Remaining, in rough order of value:

- `FE_Solution::getL2Error` is O(n²): it calls `evaluateSolution`, which linearly scans every element, once per quadrature point of every element. The element is known at the call site, so passing it in makes this linear.
- `FE_Mesh2D::evaluateSolution` returns `0.0` when no element claims the point, silently. The bounding-box reject before the affine test uses no tolerance, so a point exactly on an element boundary is at the mercy of rounding. Not observed to bite — a test over every vertex at `p = 1, 2, 3` passes — but it is luck rather than design.
- The semilinear path has an unresolved `// TODO` at `source/FE_Solution.cpp:124` for applying BCs inside the Picard loop, and no test beyond "it converges to something finite". No manufactured solution exists to check it against.
- `source/Element.cpp` is an empty file.
- Style is still inconsistent (see above); a `.clang-format` would settle it.
- `FEM.out`, a 4.6 MB Windows binary, is still in git history and so still in the clone size. Removing it needs a history rewrite and a force push, which was considered and deliberately deferred.

Extensions the dissertation itself proposes (§8), in its own order of preference: Neumann and mixed boundary conditions; hp-adaptivity driven by a posteriori error estimates; fully nonlinear problems.
