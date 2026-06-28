# Native Windows build of the cuPythia kernels (MSVC + CUDA + CMake).
# Usage:  powershell -ExecutionPolicy Bypass -File build.ps1 [-Arch <archs>] [-CudaPath <dir>]
#   -Arch     GPU architecture(s) for CMAKE_CUDA_ARCHITECTURES. Default "native" (auto-detect this
#             machine's GPU, Pascal..Blackwell). Examples: 60 (Pascal P100), 61 (GTX 10xx),
#             "all-major" (portable fatbinary), "60;70;80;90;120" (explicit list).
#   -CudaPath Specific CUDA toolkit root. Default = newest installed. For Pascal/Volta use a
#             CUDA <= 12.9 toolkit (CUDA 13 removed sm_60/61/70). "native" needs a GPU at configure
#             time; on a headless host pass an explicit -Arch. See PORTABILITY.md.
# Auto-locates Visual Studio Build Tools and its bundled CMake.
param(
  [string]$Arch = "native",
  [string]$CudaPath = ""
)
$ErrorActionPreference = "Stop"
$here = if ($PSScriptRoot) { $PSScriptRoot } else { (Get-Location).Path }

$vsw = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
$vs  = & $vsw -latest -products * -property installationPath
if (-not $vs) { Write-Error "Visual Studio Build Tools not found (install the C++ workload)."; exit 1 }
$vcvars = Join-Path $vs "VC\Auxiliary\Build\vcvars64.bat"
$cmake  = Join-Path $vs "Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
if (-not (Test-Path $cmake)) { $cmake = "cmake" }   # fall back to a cmake on PATH

if ($CudaPath) {
  $nvcc = Join-Path $CudaPath "bin\nvcc.exe"
  if (-not (Test-Path $nvcc)) { Write-Error "nvcc not found at $nvcc (check -CudaPath)."; exit 1 }
} else {
  $cudaRoot = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA"
  $cuda = Get-ChildItem $cudaRoot -ErrorAction SilentlyContinue | Sort-Object Name -Descending | Select-Object -First 1
  if (-not $cuda) { Write-Error "CUDA toolkit not found under $cudaRoot (use -CudaPath)."; exit 1 }
  $nvcc = Join-Path $cuda.FullName "bin\nvcc.exe"
}

Write-Host "VS:    $vs"
Write-Host "CMake: $cmake"
Write-Host "CUDA:  $nvcc"
Write-Host "Arch:  $Arch"

$build = Join-Path $here "build-win"
# Run configure+build inside the MSVC environment so cl.exe + nmake are found.
$cmd = "`"$vcvars`" >nul 2>&1 && set `"CUDACXX=$nvcc`" && " +
       "`"$cmake`" -S `"$here`" -B `"$build`" -G `"NMake Makefiles`" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=`"$Arch`" && " +
       "`"$cmake`" --build `"$build`""
cmd /c $cmd
if ($LASTEXITCODE -ne 0) { Write-Error "build failed ($LASTEXITCODE)"; exit 1 }
Write-Host "`nOK -> executables in $build (e.g. .\build-win\mc_pi.exe)"
