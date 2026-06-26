# Native Windows build of the cuPythia kernels (MSVC + CUDA + CMake).
# Usage:  powershell -ExecutionPolicy Bypass -File build.ps1
# Auto-locates Visual Studio Build Tools, its bundled CMake, and the newest CUDA.
$ErrorActionPreference = "Stop"
$here = if ($PSScriptRoot) { $PSScriptRoot } else { (Get-Location).Path }

$vsw = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
$vs  = & $vsw -latest -products * -property installationPath
if (-not $vs) { Write-Error "Visual Studio Build Tools not found (install the C++ workload)."; exit 1 }
$vcvars = Join-Path $vs "VC\Auxiliary\Build\vcvars64.bat"
$cmake  = Join-Path $vs "Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
if (-not (Test-Path $cmake)) { $cmake = "cmake" }   # fall back to a cmake on PATH

$cudaRoot = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA"
$cuda = Get-ChildItem $cudaRoot -ErrorAction SilentlyContinue | Sort-Object Name -Descending | Select-Object -First 1
if (-not $cuda) { Write-Error "CUDA toolkit not found under $cudaRoot."; exit 1 }
$nvcc = Join-Path $cuda.FullName "bin\nvcc.exe"

Write-Host "VS:    $vs"
Write-Host "CMake: $cmake"
Write-Host "CUDA:  $nvcc"

$build = Join-Path $here "build-win"
# Run configure+build inside the MSVC environment so cl.exe + nmake are found.
$cmd = "`"$vcvars`" >nul 2>&1 && set `"CUDACXX=$nvcc`" && " +
       "`"$cmake`" -S `"$here`" -B `"$build`" -G `"NMake Makefiles`" -DCMAKE_BUILD_TYPE=Release && " +
       "`"$cmake`" --build `"$build`""
cmd /c $cmd
if ($LASTEXITCODE -ne 0) { Write-Error "build failed ($LASTEXITCODE)"; exit 1 }
Write-Host "`nOK -> executables in $build (e.g. .\build-win\mc_pi.exe)"
