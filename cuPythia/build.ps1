# Native Windows build of the cuPythia kernels (MSVC + CUDA + CMake).
# Zero-config by default: auto-detects this machine's GPU (nvidia-smi) and picks a CUDA toolkit that
# can compile for it (so a Pascal box with both CUDA 13 and 12.9 automatically uses 12.9).
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

# Auto-detect EVERY GPU's arch (when -Arch is the default "native") so we can pick a CUDA toolkit that
# compiles for all of them and build a fatbinary if they differ -- e.g. an A100 + RTX 4090 -> {80,89},
# or a Pascal box with both CUDA 13 and 12.9 -> picks 12.9.
function Get-GpuArchs {
  $smi = (Get-Command nvidia-smi -ErrorAction SilentlyContinue).Source
  if (-not $smi -and (Test-Path "C:\Windows\System32\nvidia-smi.exe")) { $smi = "C:\Windows\System32\nvidia-smi.exe" }
  if (-not $smi) { return @() }
  $caps = & $smi --query-gpu=compute_cap --format=csv,noheader 2>$null
  if (-not $caps) { return @() }
  return @($caps | ForEach-Object { ($_ -replace '[ .]','') } | Where-Object { $_ } | Sort-Object -Unique)
}
$detected = if ($Arch -eq "native") { Get-GpuArchs } else { @() }
if ($detected.Count) { Write-Host "Detected GPU arch(s) -> $(( $detected | ForEach-Object { "sm_$_" }) -join ', ')" }

if ($CudaPath) {
  $nvcc = Join-Path $CudaPath "bin\nvcc.exe"
  if (-not (Test-Path $nvcc)) { Write-Error "nvcc not found at $nvcc (check -CudaPath)."; exit 1 }
} else {
  $cudaRoot = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA"
  # VERSION-aware sort (newest first) -- a plain string sort puts v9.0 above v12.9.
  $cudas = Get-ChildItem $cudaRoot -Directory -ErrorAction SilentlyContinue |
    Sort-Object @{ Expression = { try { [version]($_.Name -replace '^[vV]','') } catch { [version]'0.0' } } } -Descending
  if (-not $cudas) { Write-Error "CUDA toolkit not found under $cudaRoot (use -CudaPath)."; exit 1 }
  if ($detected.Count) {
    # newest toolkit that can EMIT EVERY detected arch (skips CUDA 13.x for any Pascal/Volta in the set)
    $nvcc = $null
    foreach ($c in $cudas) {
      $n = Join-Path $c.FullName "bin\nvcc.exe"; if (-not (Test-Path $n)) { continue }
      $have = & $n --list-gpu-arch 2>$null
      $ok = $true; foreach ($a in $detected) { if ($have -notcontains "compute_$a") { $ok = $false; break } }
      if ($ok) { $nvcc = $n; break }
    }
    if (-not $nvcc) {
      $hint = if (($detected | Where-Object { "60","61","70" -contains $_ })) { " Pascal/Volta need CUDA <= 12.9 (CUDA 13 removed them) -- install one or pass -CudaPath." } else { "" }
      Write-Error "No installed CUDA toolkit can compile for arch set { $($detected -join ',') }.$hint See PORTABILITY.md."; exit 1
    }
    $Arch = $detected -join ';'   # one arch, or a fatbinary list for a mixed-GPU box
  } elseif ($Arch -eq "native") {
    Write-Error "Could not detect a GPU (nvidia-smi). For a headless/CI build pass an explicit -Arch (e.g. -Arch 80, or -Arch all-major for a portable fatbinary). See PORTABILITY.md."; exit 1
  } else {
    $nvcc = Join-Path $cudas[0].FullName "bin\nvcc.exe"   # explicit -Arch given -> newest toolkit
  }
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
