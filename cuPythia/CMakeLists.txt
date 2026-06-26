# Cross-platform build for the cuPythia GPU kernels — identical on Windows and
# Linux. Requires CMake >= 3.20, a CUDA toolkit, and a host compiler (MSVC on
# Windows, g++/clang on Linux).
#
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build
#
# Override the GPU target with -DCMAKE_CUDA_ARCHITECTURES=90 (Hopper), 80, etc.
cmake_minimum_required(VERSION 3.20)

# VS 2022/2026 (MSVC 14.4x/14.5x) can be newer than the CUDA host-compiler
# allow-list. The flag must be present during CMake's CUDA compiler check, i.e.
# BEFORE project(), so set it on the CUDA flags here.
if(CMAKE_HOST_WIN32 AND NOT DEFINED CMAKE_CUDA_FLAGS)
  set(CMAKE_CUDA_FLAGS "-allow-unsupported-compiler" CACHE STRING "" FORCE)
endif()

project(cuPythia LANGUAGES CXX CUDA)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 120)   # Blackwell / RTX 5050 by default
endif()
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 14)
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
  11_vegas_importance/vegas.cu)

foreach(src ${CUPYTHIA_KERNELS})
  get_filename_component(name ${src} NAME_WE)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE Threads::Threads)
endforeach()
