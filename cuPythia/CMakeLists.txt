# Cross-platform build for the cuPythia GPU kernels — identical on Windows and
# Linux. Requires CMake >= 3.20, a CUDA toolkit, and a host compiler (MSVC on
# Windows, g++/clang on Linux).
#
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build
#
# GPU architecture (portable: Pascal sm_60 through Blackwell sm_120):
#   default              -> "native" (needs CMake >= 3.24 AND a GPU visible at configure time): auto-
#                           detects this machine's GPU, so a Pascal box on CUDA <= 12.9 builds sm_60/61
#                           with no extra flags. HEADLESS / CI / no-GPU build hosts MUST pass an explicit
#                           -DCMAKE_CUDA_ARCHITECTURES (native cannot detect a GPU and will error).
#   -DCMAKE_CUDA_ARCHITECTURES=60        Pascal P100 (or 61 for GTX 10-series)
#   -DCMAKE_CUDA_ARCHITECTURES=all-major portable fatbinary, every major arch the TOOLKIT supports (CMake >= 3.23)
#   -DCMAKE_CUDA_ARCHITECTURES="60;70;80;90;120"   explicit list (works on any CMake >= 3.18)
# NOTE: CUDA 13.x removed Pascal/Volta; to target sm_60/61/70 use CUDA <= 12.9 (see PORTABILITY.md).
cmake_minimum_required(VERSION 3.20)

# VS 2022/2026 (MSVC 14.4x/14.5x) can be newer than the CUDA host-compiler
# allow-list. The flag must be present during CMake's CUDA compiler check, i.e.
# BEFORE project(), so set it on the CUDA flags here.
if(CMAKE_HOST_WIN32 AND NOT DEFINED CMAKE_CUDA_FLAGS)
  set(CMAKE_CUDA_FLAGS "-allow-unsupported-compiler" CACHE STRING "" FORCE)
endif()

project(cuPythia LANGUAGES CXX CUDA)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
    set(CMAKE_CUDA_ARCHITECTURES native)  # auto-detect the local GPU (Pascal..Blackwell)
  else()
    set(CMAKE_CUDA_ARCHITECTURES 120)     # CMake < 3.24 has no "native": pick your arch with -D
  endif()
endif()
set(CMAKE_CXX_STANDARD 17)    # CUB (kernel 10) requires C++17 -> CUDA toolkit >= 11.0
set(CMAKE_CUDA_STANDARD 17)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Threads REQUIRED)   # std::thread in the multi-GPU kernel

set(CUPYTHIA_KERNELS
  00_toolchain_check/mc_pi.cu
  01_xsec_ee_mumu/xsec_ee_mumu.cu
  02_qcd_2to2_me/qcd_2to2.cu
  03_fused_resident_me/fused_gg2gg.cu
  04_multi_gpu/multi_gpu_gg2gg.cu
  05_reproducible_events/repro_events.cu
  06_fp32_mixed/fp32_gg2gg.cu
  07_unweighting_lhe/unweight_lhe.cu
  08_qcd_process_library/qcd_library.cu
  09_neutrino_dis/nu_dis.cu
  10_cub_compaction/cub_compaction.cu
  11_vegas_importance/vegas.cu
  12_phase_space/phase_space.cu
  13_rescattering/rescatter.cu
  14_parton_shower/shower.cu)

foreach(src ${CUPYTHIA_KERNELS})
  get_filename_component(name ${src} NAME_WE)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE Threads::Threads)
endforeach()
