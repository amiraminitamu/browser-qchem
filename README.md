# browser-qchem

Electronic structure calculations inside a web browser, using compiled
scientific libraries and no remote calculation API.

## First milestone

Compile **libcint v6.1.2** to WebAssembly and evaluate the overlap between two
normalized hydrogen STO-3G 1s orbitals. The page accepts an H-H separation in
bohr and compares libcint's result against an independent analytical formula.
This is an integral-engine demonstration, not an HF calculation yet.

The upstream library is a git submodule pinned to commit
`8d13863ff481cea27efea5e56c9e4d352cdb8f80`; its sources are unmodified.
The initial build uses scalar C and one thread. There is no Python runtime
inside the browser. The small calculation runs on the browser's main thread
for this first demo; move the engine into a worker before adding SCF workloads.

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
python3 -m http.server 8000 --bind 127.0.0.1 --directory build/wasm/dist
```

Open **http://localhost:8000**. At 1.4 bohr, the overlap should be about
**0.6593182061** and the page should report **PASS**. Changing the separation
recalculates the overlap by calling the compiled libcint code.

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

Initial validation passed on Linux x86_64 with Emscripten 4.0.10, Node 24,
and headless Chromium 134. The native comparison used PySCF 2.14.0.
The macOS setup commands are provided but have not been executed on a Mac here.

An optional native/PySCF comparison checks the wrapper's basis and normalization
conventions against PySCF:

```bash
python3 -m pip install pyscf
cmake -S . -B build/native -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build/native --parallel 4
python3 tests/check_pyscf.py build/native/libqc_core.dylib
```

On Linux the library is `build/native/libqc_core.so`. PySCF is a development
reference only. Agreement with PySCF checks packing and normalization, not
independence of the integral algorithm, because PySCF also uses libcint.
The analytical formula provides the independent check for this s-only case.
Higher angular momenta and other integral types are not validated by this demo.

## Files

| Path | Purpose |
| --- | --- |
| `vendor/libcint` | Pinned upstream integral library |
| `src/overlap.c` | Build libcint inputs and expose one overlap calculation |
| `CMakeLists.txt` | Build native or WebAssembly targets |
| `scripts/build_wasm.sh` | Configure and compile the browser module |
| `web/` | Small local test interface and analytical reference |
| `tests/` | WebAssembly and optional native/PySCF numerical checks |

## Next milestone

Expose general molecular and basis data, then validate overlap, kinetic energy,
nuclear attraction, and electron repulsion integrals for water before adding RHF.

## License and references

The project license is retained in `LICENSE`. The libcint dependency retains its
own license in `vendor/libcint/LICENSE`.

- [libcint source and documentation](https://github.com/sunqm/libcint/tree/v6.1.2)
- Q. Sun, *J. Comput. Chem.* **36**, 1664-1671 (2015),
  [doi:10.1002/jcc.23981](https://doi.org/10.1002/jcc.23981).
- W. J. Hehre, R. F. Stewart, J. A. Pople, *J. Chem. Phys.* **51**, 2657-2664
  (1969), [doi:10.1063/1.1672392](https://doi.org/10.1063/1.1672392).
- [Emscripten installation](https://emscripten.org/docs/getting_started/downloads.html)
