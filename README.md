# browser-qchem

Electronic structure calculations inside a web browser, using compiled
scientific libraries and no remote calculation API.

## Molecular integral milestone

Compute the complete overlap S, kinetic T, nuclear attraction V, and electron
repulsion tensor for editable molecular geometries. The default example is
water/STO-3G. Inspect S, T, V, T+V, individual ERIs, and download all numerical
results as JSON. This is an integral-engine demonstration, not an HF calculation.

The browser includes STO-3G basis parameters for H through Ne. Input may be a
standard XYZ file's contents (atom count, comment line, then coordinates) or
plain `Element x y z` lines. Select angstrom or bohr explicitly. The backend
uses bohr. Electron count and spin are not required for these integrals.

The upstream library is a git submodule pinned to commit
`8d13863ff481cea27efea5e56c9e4d352cdb8f80`; its sources are unmodified.
The build uses scalar C and one computation thread inside a browser worker.
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
python3 -m http.server 8000 --bind 127.0.0.1 --directory build/wasm/dist
```

Open **http://localhost:8000**. The default water calculation should show:

- 5 shells and 7 AOs.
- Three 7 x 7 one-electron matrices and 2401 electron repulsion integrals.
- Nuclear repulsion energy **9.188258417746 Eh**.
- **PASS: normalization and symmetry checks.**

The live PASS message checks normalization and integral symmetries. It is not
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

To compare all integrals directly with an installed PySCF, including H through Ne
atoms, H2, water, skew water, methane, ammonia, HF, LiH, and BeH2:

```bash
python3 -m pip install pyscf
cmake -S . -B build/native -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build/native --parallel 4
python3 tests/check_molecular.py --native build/native/libqc_core.dylib
python3 tests/check_molecular.py --wasm
```

On Linux the native library is `build/native/libqc_core.so`. The native test
also exercises custom s/p/d shells. PySCF is a development reference only.
Agreement checks basis normalization, shell ordering, units, signs, array layout,
and integration through our wrapper. PySCF also uses libcint, so this is not an
independent validation of libcint's integral algorithm. The original H2
analytical formula remains an independent check of that overlap case.

`python3 tests/generate_water_reference.py` regenerates the committed water
reference using PySCF directly. It never calls our wrapper.

For the optional end-to-end browser check:

```bash
python3 -m pip install playwright
python3 -m playwright install chromium
python3 tests/check_browser.py
```

This starts a temporary localhost server, loads the real WASM module in a worker,
compares every downloaded water integral against the reference, and checks matrix
and ERI inspection, bohr input, invalid geometry, recovery, and narrow layout.

Initial native and WASM validation passed on Linux x86_64 with Emscripten 4.0.10,
Node 24, PySCF 2.14.0, and headless Chromium 134. Across the 18 atom/molecule cases,
the maximum absolute integral difference was **2.85e-14**. The new molecular
milestone has not yet been executed on macOS here; the first H2 browser milestone
was confirmed on the user's Mac.

## Conventions and limits

AO order is atom order, then the listed shells. For the included STO-3G data it
is 1s for H/He and 1s, 2s, 2px, 2py, 2pz for Li through Ne. Both primitives and
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
basis catalog currently exposes only STO-3G s and p functions. `qc_compute`
returns an opaque result handle; read it with `qc_nao`, `qc_enuc`, `qc_data`,
then release it with `qc_free`. Inputs remain owned by the caller. No function
retains an input buffer. A null result is accompanied by a `qc_error` message.

## Files

| Path | Purpose |
| --- | --- |
| `vendor/libcint` | Pinned upstream integral library |
| `src/overlap.c` | Original H2 overlap regression path |
| `src/integrals.c`, `src/integrals.h` | General shell inputs and complete integral evaluation |
| `CMakeLists.txt` | Build native or WebAssembly targets |
| `scripts/build_wasm.sh` | Configure and compile the browser module |
| `web/basis.mjs`, `web/molecule.mjs` | STO-3G catalog, input parsing, and AO labels |
| `web/engine.mjs`, `web/worker.mjs` | WASM memory interface and background computation |
| `web/index.html`, `web/app.mjs` | Matrix inspection, ERI selection, and JSON export |
| `tests/` | Numerical comparisons, PySCF reference, and browser checks |

## Next milestone

Build RHF on the validated S, T, V, and ERIs, with orthogonalization, SCF iteration,
convergence reporting, and total energies validated against reference calculations.

## License and references

The project license is retained in `LICENSE`. The libcint dependency retains its
own license in `vendor/libcint/LICENSE`.

- [libcint source and documentation](https://github.com/sunqm/libcint/tree/v6.1.2)
- [PySCF STO-3G source](https://github.com/pyscf/pyscf/blob/v2.14.0/pyscf/gto/basis/sto-3g.dat).
  The included numerical basis parameters are transcribed from PySCF 2.14.0;
  no Python is needed to use them in the browser.
- Q. Sun, *J. Comput. Chem.* **36**, 1664-1671 (2015),
  [doi:10.1002/jcc.23981](https://doi.org/10.1002/jcc.23981).
- W. J. Hehre, R. F. Stewart, J. A. Pople, *J. Chem. Phys.* **51**, 2657-2664
  (1969), [doi:10.1063/1.1672392](https://doi.org/10.1063/1.1672392).
- [Emscripten installation](https://emscripten.org/docs/getting_started/downloads.html)
