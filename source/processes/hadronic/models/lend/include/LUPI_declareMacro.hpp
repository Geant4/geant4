/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef LUPI_declare_macro_hpp_included
#define LUPI_declare_macro_hpp_included

#include <LUPI_defines.hpp>

// Default, if LUPI_HIP_INLINE is not defined, is to use an attribute function
// to inform HIP not to inline the function.  
// This is quite useful for the publicly installed header files for code robustness
// However, when compiling the source, one should be able to disable
// this within the library itself for faster code.  To  do so
// the define -DLUPI_HIP_INLINE can be added to the compiler flags, and
// the define will evaluate to nothing, so the compiler is welcome to do 
// its optimizations.
   
#ifdef LUPI_HIP_INLINE
    #define LUPI_HIP_INLINE_ATTRIBUTE
#else
    #define LUPI_HIP_INLINE_ATTRIBUTE  __attribute__ ((noinline))
#endif

#define gpuErrorCheck(ans) { gpuAssert((ans), __FILE__, __LINE__); }

#if defined(__HIP_DEVICE_COMPILE__) || defined(__CUDA_ARCH__)
    #define LUPI_ON_GPU 1
#endif

#ifdef __CUDACC__
    #include <cstdio>
    #define LUPI_HOST __host__
    #define LUPI_DEVICE __device__
    #define LUPI_HOST_DEVICE __host__ __device__
    #define LUPI_THROW(arg) printf("%s", arg)
    #define LUPI_WARP_SIZE 32
    #define LUPI_THREADID threadIdx.x
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUASSERT: %s File: %s line: %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#elif HAVE_OPENMP_TARGET
    #define LUPI_HOST 
    #define LUPI_DEVICE 
    #define LUPI_HOST_DEVICE 
    #define LUPI_THROW(arg) printf("%s", arg)
    #define LUPI_WARP_SIZE 1
    #define LUPI_THREADID 
inline void gpuAssert(int code, const char *file, int line, bool abort=true) {}
#elif defined(__HIP__)
    #include <hip/hip_version.h>
    #include <hip/hip_runtime.h>
    #include <hip/hip_runtime_api.h>
    #include <hip/hip_common.h>

    #define LUPI_HOST __host__
    #define LUPI_DEVICE __device__
    #define LUPI_HOST_DEVICE LUPI_HIP_INLINE_ATTRIBUTE __host__ __device__
    #define LUPI_THROW(arg) 
    #define LUPI_WARP_SIZE 1
    #define LUPI_THREADID hipThreadIdx_x
inline void gpuAssert(hipError_t code, const char *file, int line, bool do_abort=true)
{
    if (code == hipSuccess) { return; }
    printf("GPUassert code %d: %s %s %d\n", code, hipGetErrorString(code), file, line);
    if (do_abort) { abort(); }
}

#else
    #define LUPI_HOST
    #define LUPI_DEVICE 
    #define LUPI_HOST_DEVICE
    #define LUPI_THROW(arg) throw arg
    #define LUPI_WARP_SIZE 1
    #define LUPI_THREADID 
inline void gpuAssert(LUPI_maybeUnused int code, LUPI_maybeUnused const char *file, LUPI_maybeUnused int line, LUPI_maybeUnused bool abort=true) {}
#endif

#endif      // End of LUPI_declare_macro_hpp_included
