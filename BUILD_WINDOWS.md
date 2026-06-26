# Building cuPythia on Windows 11 (native, no WSL)

Stock Pythia 8 has **no native Windows build** — its `configure` is a Unix bash
script and the code assumes a POSIX toolchain. cuPythia's GPU kernels, by
contrast, build and run **natively on Windows 11** with the *same source* and the
*same validations* as Linux. The only portability shim required was an `M_PI`
fallback (MSVC's `<cmath>` omits `M_PI` unless `_USE_MATH_DEFINES` is set).

**Verified:** MSVC 14.50 (VS Build Tools 2026) + CUDA 13.1 + RTX 5050 — all 9
kernels configure, build, and `VALIDATION: PASS` natively.

## Requirements
- **Visual Studio Build Tools** (the "Desktop development with C++" workload) —
  provides `cl.exe` and a bundled CMake.
- **NVIDIA CUDA Toolkit** (>= 12.8 for Blackwell / `sm_120`).
- An NVIDIA GPU + driver.

## One command
```powershell
cd cuPythia\cuPythia
powershell -ExecutionPolicy Bypass -File build.ps1
.\build-win\mc_pi.exe
.\build-win\qcd_library.exe
```
`build.ps1` auto-locates Build Tools, its CMake, and the newest CUDA toolkit, then
configures and builds (NMake) inside the MSVC environment.

## Or CMake directly
```powershell
# from an "x64 Native Tools Command Prompt for VS":
set "CUDACXX=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.1\bin\nvcc.exe"
cmake -S cuPythia\cuPythia -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```
Override the GPU target with `-DCMAKE_CUDA_ARCHITECTURES=90` (Hopper), `80`, etc.
The same `CMakeLists.txt` builds on Linux (g++/clang + nvcc) unchanged.

## Notes
- **Clone to a local path** (e.g. `C:\src\cuPythia`). `cmd.exe` cannot use a
  `\\wsl$\...` UNC path as a working directory; `build.ps1` works around it with
  absolute paths, but a local clone is cleanest.
- VS 2026's MSVC (14.5x) is newer than CUDA's host-compiler allow-list, so the
  build passes `-allow-unsupported-compiler` automatically.
- **Pythia itself** (`pythia8317/`) still builds via its Unix toolchain — use WSL
  for the full library. The GPU **kernels** are what run natively on Windows; a
  native MSVC/CMake port of all of Pythia 8 is a larger, separate effort.
