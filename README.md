# Lorentzian

<img src="docs/lorentzian-icon.jpg" width="88" alt="Lorentzian mark: a copper Lorentzian peak on navy">

Electronic structure calculations inside a web browser, using compiled
scientific libraries and no remote calculation API. The name is the usual
English adjective from Lorentz — the same pattern as Gaussian, and the competing
spectroscopic lineshape.

The packaged workbench is at
**https://amiraminitamu.github.io/browser-qchem/**.
GitHub Pages serves the WASM build from the `gh-pages` branch. Enable it once
under **Settings → Pages → Deploy from a branch → `gh-pages` / `/ (root)`**.
Later pushes to `main` rebuild that branch through `.github/workflows/pages.yml`.

## Molecular workspace

Run RHF, unrestricted Hartree-Fock (UHF), and small conventional RHF-MP2 calculations directly in the browser. Enter a
geometry and molecular charge, follow live SCF iterations, inspect orbital
energies, density and Fock matrices, and export numerical results as JSON. UHF adds
spin multiplicity, separate alpha/beta orbitals and matrices, and spin diagnostics.
MP2 adds the all-electron correlation and total energies on an RHF reference.
Choose STO-3G, 6-31G, or cc-pVDZ for H through Ne. The workspace includes a
rotatable 3D molecule, XYZ upload, atom picking, live SCF and PES plots, matrix
inspection, and numerical exports. The default water/STO-3G calculation converges
in 9 iterations to
**-74.963063129729 Eh** with DIIS enabled.

The solver uses libcint integrals and Eigen diagonalization and linear solves.
The default convergence criteria are energy change < 1e-10 Eh, AO density RMS
change < 1e-8, and orthonormal-basis commutator RMS < 1e-8. Results that reach
the iteration limit are labeled **NOT CONVERGED**, including in JSON metadata.

## Molecular integral foundation

Compute the complete overlap S, kinetic T, nuclear attraction V, and electron
repulsion tensor for editable molecular geometries. The default example is
water/STO-3G. Inspect S, T, V, T+V, individual ERIs, and download all numerical
results as JSON. Select **Integrals only** to use the earlier integral calculator.

The browser bundles all three basis catalogs locally, including general
contractions and spherical d functions in cc-pVDZ. Input may be a
standard XYZ file's contents (atom count, comment line, then coordinates) or
plain `Element x y z` lines. Select angstrom or bohr explicitly. The backend
uses bohr. RHF requires a positive even electron count and assumes spin 0;
electron count and spin are not required for integrals-only calculations.

libcint is pinned to commit `8d13863ff481cea27efea5e56c9e4d352cdb8f80` (v6.1.2).
Eigen is pinned to `3147391d946bb4b6c68edd901f2add6ac1f31f8c` (3.4.0).
Both are unmodified git submodules. Updating from the integral milestone requires
`git submodule update --init --recursive` to download Eigen's headers.
The build uses C/C++ and one computation thread inside a browser worker.
There is no Python runtime or remote calculation API in the browser. No special
cross-origin isolation headers are needed for this single-worker implementation.

## Set up on macOS

Git, Python 3, CMake >= 3.20, a Make-compatible build tool, and Emscripten
are needed to build. Git already works on your machine. If needed, install
Apple's command-line tools with `xcode-select --install`.

With conda, create a dedicated build environment:

```bash
conda create -n browser-qchem -c conda-forge python=3.12 cmake -y
conda activate browser-qchem
```

From the repository root, install Emscripten **4.0.10** into an ignored tools
directory. No administrator access or Homebrew is required:

```bash
mkdir -p .tools
git clone https://github.com/emscripten-core/emsdk.git .tools/emsdk
./.tools/emsdk/emsdk install 4.0.10
./.tools/emsdk/emsdk activate 4.0.10
source .tools/emsdk/emsdk_env.sh
emcc --version
```

Downloads are substantial on the first install. If an SDK is already installed,
activate it and source its `emsdk_env.sh` instead. Version 4.0.10 is the initial
tested toolchain rather than a requirement to track the latest SDK.

## Build and open the page

From the repository root:

```bash
git submodule update --init --recursive
bash scripts/build_wasm.sh
node tests/check_wasm.mjs
node tests/check_molecular.mjs
node tests/check_rhf.mjs
node tests/check_workspace.mjs
node tests/check_methods.mjs
node tests/check_state.mjs
python3 scripts/serve.py --directory build/wasm/dist
```

Open **http://localhost:8000**. The default water calculation should show:

- 5 shells and 7 AOs.
- Three 7 x 7 one-electron matrices and 2401 electron repulsion integrals.
- Nuclear repulsion energy **9.188258417746 Eh**.
- Total RHF energy **-74.963063129729 Eh**.
- **CONVERGED: RHF/STO-3G.**

In integrals-only mode, the live PASS message checks normalization and integral symmetries. It is not
a runtime PySCF comparison for arbitrary entered geometries. The included
numerical tests compare full arrays against PySCF. On updating from the first
demo, rebuild and reload the browser with Cmd+Shift+R to refresh the JS and WASM.

The local Python server only serves static files; all integral evaluation
happens in the browser. No external assets or calculation services are used.
Stop the server with Ctrl+C. Use HTTP rather than opening the HTML as a file.

For a new terminal session, reactivate the tools before rebuilding:

```bash
conda activate browser-qchem
source .tools/emsdk/emsdk_env.sh
```

Re-run the build after changing the C code or webpage. CMake copies the webpage
into `build/wasm/dist`; do not edit the generated copies there.

## Settings and build consistency

Method changes commit immediately on the select's change event. Presets recommend
a method only until the user explicitly selects one. They still load the
appropriate charge and multiplicity. **Reset** restores all defaults, including
RHF and the water geometry.

Each job captures one immutable request before controls are disabled. The same
snapshot supplies validation, labels, scan metadata and the worker request.
An old scan cannot rewrite a different running calculation. Inspecting a saved
point preserves newer method, basis, charge and spin edits, with an explicit
saved-result label. Edits invalidate
late replies; wrong-method, wrong-basis or wrong-charge results are rejected.
The worker requires an explicit method and a matching protocol/build handshake.

The complete editable draft is saved locally in the browser under
`browser-qchem.settings.v2`. Reloading or reopening restores geometry, method,
basis, task, charge, spin and scan settings together. Only settings are saved;
numerical result arrays are not stored. Storage-disabled browsers still work.
A page first opened while the user is editing waits for Run instead of starting
an automatic calculation with partially edited settings.

Every WASM build packages the complete local module graph with one content hash.
The page, styles, worker, module imports and WASM request use that same version.
Use `python3 scripts/serve.py` to serve with no-store headers. Stop any previous
server before starting it. After this update, perform one Cmd+Shift+R refresh to
replace the previously cached, unversioned page.

Additional checks:

```bash
node tests/check_state.mjs
python3 tests/check_package_web.py
python3 tests/check_browser_state.py
```

The controller/worker regression covers explicit methods versus presets,
active-job scan inspection, stale replies, native-select event order, full
reset, coherent draft restoration, dispatch and build mismatches. The real-WASM
browser regression also covers repeated RHF/UHF/MP2 changes, reload and reopening.

## Interface and navigation

The workbench uses a canvas-first layout: a single compute pipeline above a full-height
molecular stage and a unified results dock. The app stays inside the viewport;
plots, tables, and drawers scroll within their own panes.

- The **Molecule pill** contains examples, XYZ import, and reset. **XYZ** opens
  the coordinate editor over the left edge of the canvas.
- The **Task pill** selects single-point energy, bond-length scan, or molecular
  integrals. Scan ranges and fragment controls appear inside this popover.
- **Method / basis**, **charge and spin**, and **Solver** settings share the
  pipeline. Solver controls adapt to RHF, UHF, MP2, and integral-only tasks.
- **Run** becomes **Stop** during computation. **Cmd+Enter** or **Ctrl+Enter**
  runs the selected calculation. Completed calculations appear in **Jobs**.
- **Export** is the home for geometry XYZ, result JSON, and scan CSV/XYZ/JSON.
  Canvas PNG export sits with labels, fit, and expand in the floating viewer toolbar.
- Drag the divider to resize the stage from 35% to 70%, or focus the divider
  and use arrow keys, Home, and End. Double-click restores the 55% split.
- The results dock keeps energy, status, cycles, and frontier orbital energies
  visible above its tabs. **Conventions** includes numerical and method diagnostics.
- Tabs support arrow keys, Home, and End. SCF chart points highlight iteration
  rows. Orbital rows expose selection events for future isosurface rendering;
  the current implementation displays the selected orbital and its energy.
- **Geometry** contains Cartesian coordinates, internuclear distances, and
  **Copy XYZ**. **Esc** closes popovers, jobs, the coordinate drawer, or viewer expansion.
- Below 800 px the stage and dock stack vertically. The same pipeline controls
  remain available; no duplicate setup or run bar is added.

Task/tab registries, method-specific solver groups, canvas overlay slots, and
split/drawer behavior live in `web/workbench.mjs`. Calculation settings and
worker requests remain in `web/app.mjs` and `web/settings.mjs`.

### Ions and open-shell examples

The example selector includes systems with electronic recommended settings loaded together when no method has been explicitly selected:

| Example | Charge | Electrons | Multiplicity | Default method |
| --- | --- | --- | --- | --- |
| H2+ | +1 | 1 | 2 (doublet) | UHF |
| HeH, neutral | 0 | 3 | 2 (doublet) | UHF |
| HeH+ | +1 | 2 | 1 (singlet) | RHF |

HeH is neutral; HeH+ is the cation. Selecting either open-shell example recommends
UHF and sets the core-Hamiltonian guess. An explicitly selected method takes
precedence over subsequent example recommendations. You can edit charge and multiplicity for
other ions, including anions. The formula and electron count update with charge.
RHF and RHF-MP2 continue to require a positive even electron count, so H2+
requires UHF. Geometry-only XYZ input does not carry charge or multiplicity;
set those fields explicitly after importing a file.

The ion regression checks include H2+, neutral HeH, HeH+ with UHF, and H- with
UHF in STO-3G, 6-31G and cc-pVDZ, plus all-electron RHF-MP2 for HeH+.
The full reference suite covers **56 single points and 27 ion/open-shell PES
points** against PySCF. Scan CSV and XYZ now include charge and electron count
in addition to method, multiplicity and energy. Browser tests verify preset
settings, superscript charges, charged scans and exports, invalid spin states,
method rejection for odd-electron MP2, and mobile recovery.

## Explore a molecule and scan a bond

1. Select an example or upload an XYZ file. Choose the coordinate unit and basis.
   The viewer and AO count update before you calculate. Water uses 7 AOs in
   STO-3G, 13 in 6-31G, and 24 in cc-pVDZ. Calculations above 40 AOs are disabled.
2. Drag to rotate, scroll to zoom, and use **Fit** to frame the structure.
   Atom labels match the one-based indices in the geometry input. **PNG** exports
   the current view; **Export XYZ** writes an XYZ file in angstrom.
3. Select **Bond-length scan**. Click atom A and then atom B, or type their indices.
   A stays fixed. B and every listed moving atom translate along the original AB
   direction. **Select connected fragment** proposes the fragment on B's side after
   removing AB from an inferred bond graph. Check the highlighted selection;
   cyclic structures require an explicit moving group. Inferred bonds affect
   display/selection only and do not enter the electronic structure calculation.
4. Set the start/end distance, angstrom or bohr, and 2 to 101 points. Descending
   ranges are supported. Run to calculate points sequentially in the worker.
   Every geometry is generated from the original input, avoiding accumulated drift.
5. Select a curve point or table row to inspect its geometry and SCF details.
   The viewer, editor, energy and matrices follow that point. Relative energies
   use the lowest converged point in the current scan. Unconverged points are
   marked and are not connected by the curve; errors retain an explanatory row.
6. **Stop calculation** cancels between integral/SCF steps. A synchronous WASM
   integral call finishes before cancellation can be handled. Completed points
   remain available. Export energies as CSV, all completed geometries as a
   multi-frame XYZ, or compact results and metadata as JSON. Failed geometries
   remain in XYZ/JSON with their status; CSV leaves unavailable energies blank.

These are **unrelaxed bond scans**: the other coordinates are not optimized.
Each point uses an independent guess, with the selected orbital-mixing policy for UHF. There is no state
following, orbital stability test, or guarantee of the lowest SCF solution.
Core guesses can converge to excited SCF roots, including for N2. Inspect the
curve for discontinuities and verify the electronic state before interpreting a
PES. Restricted HF also cannot describe general bond dissociation correctly.

Single-point JSON contains all integrals. Scan JSON omits the large ERI tensor at
each point, retaining one-electron matrices and electronic results. To inspect a scan
geometry's ERIs, select the point and run **Single-point energy**. Exported JSON
molecular coordinates are in bohr; XYZ coordinates are in angstrom. Units and
status accompany all scan energies. Selecting an older completed point during a
scan keeps that point selected while later points continue to calculate.

The viewer uses the bundled 3Dmol.js 2.5.5 build. Without WebGL, calculations and
atom-index scan controls still work. No external asset is required at runtime.

## Numerical checks

`node tests/check_wasm.mjs` runs the actual compiled WebAssembly module at
eight separations, including zero separation as a normalization check, and
checks invalid inputs. Its reference is the closed-form contracted s-Gaussian
overlap with the same basis coefficients, normalized as a contraction.
The absolute comparison tolerance is `1e-12`.

`node tests/check_molecular.mjs` checks every entry of the water integral arrays
against the committed PySCF 2.14.0 reference in `tests/fixtures/water-sto3g.json`.
It also checks unit conversion, H2 regression, rejected inputs, recovery, and
correct heap access after WebAssembly memory growth. No PySCF installation is
needed to run the Node checks.

`node tests/check_rhf.mjs` compares the computed water energy, density, Fock, and
orbital energies with a committed PySCF reference, then checks convergence
thresholds, DIIS off, a forced iteration limit, molecular charge, invalid input,
near-linear-dependence rejection, and recovery. All four Node checks run the
actual compiled WASM engine, not a JavaScript numerical substitute.

`node tests/check_workspace.mjs` compares water and five-point H2 scans in each
basis against committed PySCF references. It checks rigid fragment geometry,
reverse ranges, bohr units, ring selection, AO limits, cancellation during SCF,
partial exports, unconverged/failed points, and recovery. Regenerate the small
reference with `python3 tests/generate_workspace_reference.py`. Regenerate the
basis catalog with `python3 scripts/generate_basis.py` using PySCF 2.14.0.

To compare all integrals directly with an installed PySCF, including H through Ne
atoms, H2, water, skew water, methane, ammonia, HF, LiH, and BeH2:

```bash
python3 -m pip install pyscf
cmake -S . -B build/native -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build/native --parallel 4
python3 tests/check_molecular.py --native build/native/libqc_core.dylib
python3 tests/check_molecular.py --wasm
python3 tests/check_rhf.py --native build/native/libqc_core.dylib
python3 tests/check_rhf.py --wasm
python3 tests/check_basis.py
```

The basis check compares every S, T, V and ERI entry plus RHF energies and
densities for water, skew water, ammonia and methane in all three bases. In 12
cases, maximum differences were 3.56e-14 for integrals, 3.63e-13 Eh for energies,
and 1.12e-8 for density entries. The cc-pVDZ cases exercise general contractions
and all five spherical d functions. All integral entries for the 30 H-Ne
atom/basis combinations also pass, with maximum difference 4.27e-14.

On Linux the native library is `build/native/libqc_core.so`. The native test
also exercises custom s/p/d shells. PySCF is a development reference only.
Agreement checks basis normalization, shell ordering, units, signs, array layout,
and integration through our wrapper. PySCF also uses libcint, so this is not an
independent validation of libcint's integral algorithm. The original H2
analytical formula remains an independent check of that overlap case.

`python3 tests/generate_water_reference.py` regenerates the committed water
reference using PySCF directly. It never calls our wrapper.
`python3 tests/generate_rhf_reference.py` similarly regenerates the water RHF
reference using PySCF's RHF solver.

For the optional end-to-end browser check:

```bash
python3 -m pip install playwright numpy pillow
python3 -m playwright install chromium
python3 tests/check_browser.py
python3 tests/check_browser_rhf.py
python3 tests/check_browser_workspace.py
python3 tests/check_browser_interface.py
python3 tests/check_browser_methods.py
python3 tests/check_browser_ions.py
```

This starts a temporary localhost server, loads the real WASM module in a worker,
compares every downloaded water integral against the reference, and checks matrix
and ERI inspection, bohr input, invalid geometry, recovery, and narrow layout.
The RHF browser test additionally checks energy and density in downloaded JSON,
orbital energies, density/Fock matrix inspection, live iteration messages, a
forced nonconverged calculation, odd-electron rejection, and recovery.

The workspace browser check covers actual WebGL rendering, rotation, atom
picking, connected-fragment selection, PNG export, cc-pVDZ worker results,
AO preflight, XYZ upload, PES point inspection, CSV/XYZ/JSON downloads, stopping
a scan, failed-point cleanup, nonconvergence, recovery and mobile layout.

The interface browser check targets the resizable two-zone layout, canonical actions,
orbital selection, retained results after analysis changes, convergence controls,
keyboard navigation, help, viewer expansion, mobile scans, and widths from 360 to
1512 pixels. Browser interaction helpers reveal each control's popover before editing.
The canvas-first redesign has passed state/worker checks, packaging checks, and
real-WASM numerical tests. Its browser interaction and visual checks still need to run.

Initial native and WASM validation passed on Linux x86_64 with Emscripten 4.0.10,
Node 24, PySCF 2.14.0, and headless Chromium 134. Across the 18 atom/molecule cases,
the maximum absolute integral difference was **2.85e-14**. The integral and RHF
baselines were also confirmed on the user's Mac. The new basis/viewer/scan
workspace is validated on Linux and awaits its first Mac run.

The native and WASM RHF checks compare 13 cases against PySCF RHF with a
core-Hamiltonian guess: H2, stretched H2, water, skew water, methane, ammonia,
HF, LiH, BeH2, He, Ne, HeH+, and H3+. Maximum energy difference was **2.0e-13 Eh**,
maximum density-entry difference **6.7e-9**, and orbital-energy differences were
below **1.4e-9 Eh**. Tests also verify MO orthonormality, electron count,
spin-summed density idempotency, the physical Fock matrix, the generalized
eigenvalue equation, and independent recomputation of the total energy.

## UHF and conventional RHF-MP2

Select **Electronic method** independently of **Calculation type**. All methods
support single points, unrelaxed bond scans, live progress, cancellation, and
JSON exports with explicit method metadata. CSV and XYZ scan exports identify
the method and multiplicity. MP2 scans plot the correlated total energy. If the
RHF reference fails to converge, MP2 is not computed and that scan point has no
MP2 energy. Completed scan points are retained when a job is stopped.

For UHF, enter charge and multiplicity M = 2S + 1. The solver uses
Nalpha = (Ne + M - 1)/2 and Nbeta = (Ne - M + 1)/2, with nonnegative integer
populations that fit the basis. Odd electron counts require even multiplicity;
even counts require odd multiplicity. The default guess applies opposite
45-degree HOMO/LUMO rotations to alpha/beta occupied orbitals for equal spin
populations when virtual orbitals exist. Otherwise it uses the core Hamiltonian.
Select **Core Hamiltonian** to disable mixing. For stretched singlet H2, the
mixed guess allows a broken-symmetry UHF solution while a spin-symmetric core
guess can remain on the RHF solution.

UHF uses Falpha = H + J(Dalpha + Dbeta) - K(Dalpha), and analogously for beta.
Joint DIIS uses both spin commutators. Density and residual RMS are normalized
by sqrt(2) times the AO count. Both spins must pass convergence and a final
physical-Fock canonicalization check. The interface reports alpha/beta counts,
requested multiplicity, <S^2>, target S(S+1), and spin contamination. UHF is
not generally a total-spin eigenfunction. Convergence does not guarantee the
lowest or stable SCF solution; different guesses and DIIS trajectories can
select distinct stationary roots, particularly in open-shell species.

MP2 currently means **all-electron, conventional RHF-MP2** with no frozen core,
no density fitting, and no UMP2. The calculation requires converged RHF orbitals.
Each occupied orbital is transformed in four stages, retaining only O(N^3)
additional intermediate storage. The AO ERIs remain O(N^4); the transformation
cost is O(N^5). One occupied block is processed per worker step so progress and
cancellation are available between blocks. The existing **40 AO** limit remains.

The correlation energy is the sum over all occupied i,j and virtual a,b of
(ia|jb) [2(ia|jb) - (ib|ja)] / (epsilon_i + epsilon_j - epsilon_a - epsilon_b).
Separate opposite-spin and same-spin energies are reported. Nonfinite,
nonnegative, or nearly zero denominators (absolute value <= 1e-10 Eh) are
rejected rather than shifted. No virtual orbitals gives zero correlation.
MP2 can become unreliable near bond dissociation or a small RHF orbital gap.
Its displayed orbitals, density, Fock matrices and SCF history are explicitly
those of the RHF reference, not a correlated MP2 density or orbital set.

JSON keeps the existing `rhf` object for RHF and the MP2 reference. UHF uses
`uhf` with `densityAlpha`, `densityBeta`, `fockAlpha`, `fockBeta`,
`coefficientsAlpha`, `coefficientsBeta`, `orbitalEnergiesAlpha/Beta`,
`occupationsAlpha/Beta`, total `density`, and `spinDensity`. MP2 adds `mp2`
with reference, correlation, total, opposite-spin and same-spin energies,
occupied/virtual counts, gap, and minimum absolute denominator. Failed RHF
references have `mp2.status = "not_computed"` and `mp2.energy = null`.

Validation commands (PySCF is only needed to regenerate/check references):

```bash
node tests/check_methods.mjs
python3 tests/check_methods.py
python3 tests/check_methods_native.py build/native/libqc_core.dylib
python3 tests/check_browser_methods.py
python3 tests/check_browser_ions.py
```

On Linux use `build/native/libqc_core.so`. The numerical cases cover H, He+,
Li, OH, O and O2 spin states, closed-shell water, broken-symmetry H2, and MP2
for H2, water, HF, helium and ammonia across the three supported bases.
References are generated with PySCF 2.14.0. UHF checks also verify electron
counts, density idempotency, MO orthonormality, canonical residuals and physical
Fock reconstruction. The fixture records cases where PySCF's initial trajectory
chooses a different stationary root; the reported root is independently
reconverged with PySCF and is not represented as a ground-state guarantee.

## RHF algorithm and conventions

The C++ core uses a symmetric S^(-1/2) orthogonalizer and a core-Hamiltonian guess.
Its spin-summed density has occupation 2 per occupied spatial orbital:

```text
D = 2 C_occ C_occ^T
F_ij = H_ij + sum_kl D_kl [ (ij|kl) - 0.5 (ik|jl) ]
E_total = 0.5 sum_ij D_ij (H_ij + F_ij) + E_nuc
```

DIIS stores up to 8 physical Fock matrices and their orthonormal-basis
commutators, X^T (FDS-SDF) X. A singular or ill-conditioned DIIS solve falls back
to the physical Fock matrix. Convergence requires all three criteria at once:
absolute energy change, density RMS change `||D_new-D_old||_F/N`, and residual
RMS `||X^T(FDS-SDF)X||_F/N`. A final canonicalization with the physical Fock
includes a checked pure SCF step before a calculation is labeled converged.

The smallest overlap eigenvalue must be >= 1e-8. The prototype rejects nearly
linearly dependent bases rather than silently dropping AOs. RHF pairs electrons
in spatial orbitals, so odd-electron and open-shell calculations need a later
ROHF extension; UHF is available separately. SCF convergence does not establish wavefunction stability or
that the solution is the global minimum. An iteration-limited energy is exported
as an unconverged result and is labeled that way on the page.

RHF JSON includes `rhf.energy`, `electronicEnergy`, `converged`, `status`,
`charge`, `electrons`, `spin`, `iterations`, `thresholds`, `history`, `density`,
`fock`, `coefficients`, `orbitalEnergies`, and `occupations`. All matrices use
flat row-major layout. `coefficients[i*N+p]` is MO p's coefficient on AO i.
MO energies are in Hartree; occupations are 2 or 0. The C RHF API borrows the
integral handle, so free the RHF handle before freeing its integrals.

## Conventions and limits

AO order is atom order, then the listed shells. For STO-3G it is
1s for H/He and 1s, 2s, 2px, 2py, 2pz for Li through Ne. General contractions
are split into separate shell records in PySCF AO order. Display labels s1, s2,
p1, etc. denote contraction order, not assigned atomic principal quantum numbers.
Spherical d functions follow libcint order xy, yz, z2, xz, x2-y2. Both primitives and
contractions are normalized in the C wrapper. Nuclear attraction includes its
negative sign. The angstrom-to-bohr conversion uses 0.52917721092 angstrom/bohr
to match the reference calculations.

Every exported array is flat and row-major. For N AOs:

- S, T, V each have N*N values, index `i*N+j`.
- ERI has N^4 values in chemists' notation `(ij|kl)`, index `((i*N+j)*N+k)*N+l`.
- JSON indices are zero-based; the interface's AO selectors are one-based.
- Coordinates in exported molecule data are in bohr. S is dimensionless;
  T, V, ERIs, and nuclear repulsion use atomic units.

The prototype limits input to 32 atoms and **40 AOs**, with a dense ERI tensor
of at most 20.48 MB. All shell quartets are evaluated; direct Fock construction,
screening strategies, tensor compression, and parallel integral evaluation are
later optimizations. Co-located nuclei, unsupported elements, malformed input,
and excessive sizes produce explicit errors. Editing geometry clears old results.

The generic C API in `src/integrals.h` accepts arrays of nuclear charges,
coordinates, shell records `[atom_index,l,nprim]`, exponents, and contraction
coefficients. It supports one contraction per shell with l=0..2 and positive
exponents. Split multiple contractions into separate shell records. The browser
catalog exposes STO-3G, 6-31G and cc-pVDZ through spherical d functions. `qc_compute`
returns an opaque result handle; read it with `qc_nao`, `qc_enuc`, `qc_data`,
then release it with `qc_free`. Inputs remain owned by the caller. No function
retains an input buffer. A null result is accompanied by a `qc_error` message.

## Files

| Path | Purpose |
| --- | --- |
| `vendor/libcint` | Pinned upstream integral library |
| `vendor/eigen` | Pinned header-only linear algebra library |
| `src/overlap.c` | Original H2 overlap regression path |
| `src/integrals.c`, `src/integrals.h` | General shell inputs and complete integral evaluation |
| `src/rhf.cpp`, `src/rhf.h` | RHF iteration, DIIS, diagonalization, and convergence |
| `src/uhf.cpp`, `src/uhf.h` | Spin-resolved UHF, joint DIIS and spin diagnostics |
| `src/mp2.cpp`, `src/mp2.h` | Blockwise all-electron conventional RHF-MP2 |
| `CMakeLists.txt` | Build native or WebAssembly targets |
| `scripts/build_wasm.sh` | Configure and compile the browser module |
| `web/basis.mjs`, `web/molecule.mjs` | Three basis catalogs, input parsing, and AO labels |
| `web/engine.mjs`, `web/worker.mjs` | WASM memory interface and background computation |
| `web/rhf.mjs`, `web/uhf.mjs`, `web/mp2.mjs`, `web/electronic.mjs` | Method validation, memory ownership, iterations and progress |
| `web/index.html`, `web/styles.css`, `web/app.mjs` | Workspace, results inspection and exports |
| `web/workbench.mjs` | Task/tab registries, solver groups, popovers, drawers, and split resizing |
| `web/viewer.mjs`, `web/vendor/` | Local 3Dmol molecule viewer and license |
| `web/scan.mjs`, `web/plot.mjs` | Scan geometries, runner, exports and interactive SVG charts |
| `tests/` | Numerical comparisons, PySCF reference, and browser checks |

## License and references

The 3Dmol.js license and dependency notices are retained in
`web/vendor/3Dmol-LICENSE.txt`; its pinned source and checksum are recorded in
`web/vendor/README.md`.

The project license is retained in `LICENSE`. The libcint dependency retains its
own license in `vendor/libcint/LICENSE`.
Eigen retains its license files in `vendor/eigen/`; `EIGEN_MPL2_ONLY` is set
for this build.

- [libcint source and documentation](https://github.com/sunqm/libcint/tree/v6.1.2)
- [Eigen 3.4.0 source](https://gitlab.com/libeigen/eigen/-/tree/3.4.0)
- [3Dmol.js source](https://github.com/3dmol/3Dmol.js)
- [PySCF STO-3G source](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/gto/basis/sto-3g.dat).
- [PySCF 6-31G source](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/gto/basis/pople-basis/6-31G.dat)
- [PySCF cc-pVDZ source](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/gto/basis/cc-pvdz.dat)
  The included numerical basis parameters are transcribed from PySCF 2.14.0;
  no Python is needed to use them in the browser.
- Q. Sun, *J. Comput. Chem.* **36**, 1664-1671 (2015),
  [doi:10.1002/jcc.23981](https://doi.org/10.1002/jcc.23981).
- W. J. Hehre, R. F. Stewart, J. A. Pople, *J. Chem. Phys.* **51**, 2657-2664
  (1969), [doi:10.1063/1.1672392](https://doi.org/10.1063/1.1672392).
- [Emscripten installation](https://emscripten.org/docs/getting_started/downloads.html)


