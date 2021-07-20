#ifndef GCTA_CPU_H
#define GCTA_CPU_H

#if defined(__x86_64__) || (defined(_M_X64) && !defined(_M_ARM64EC)) || defined(__amd64)
  #define GCTA_ARCH_x86_64 1
#else
  #define GCTA_ARCH_x86_64 0
#endif

#if defined(__i386__) || defined(_M_IX86) || defined(_X86_) || defined(__i386)
  #define GCTA_ARCH_i386 1
#else
  #define GCTA_ARCH_i386 0
#endif

#if GCTA_ARCH_x86_64 || GCTA_ARCH_i386
  #define GCTA_CPU_x86 1
#else
  #define GCTA_CPU_x86 0
#endif

#if defined(__arm__)
  #define GCTA_ARCH_ARM 1
#else
  #define GCTA_ARCH_ARM 0
#endif

#if defined(__aarch64__) || defined(_M_ARM64) || defined(_M_ARM64EC)
  #define GCTA_ARCH_ARM64 1
#else
  #define GCTA_ARCH_ARM64 0
#endif

#if GCTA_ARCH_ARM || GCTA_ARCH_ARM64
  #define GCTA_CPU_ARM 1
#else
  #define GCTA_CPU_ARM 0
#endif

#if GCTA_CPU_x86
  #ifndef EIGEN_USE_MKL_ALL
  #define EIGEN_USE_MKL_ALL
  #endif
  #include <mkl.h>
#else
  #ifndef EIGEN_USE_BLAS
  #define EIGEN_USE_BLAS
  #endif
  #include <cblas.h> 
  #include <lapack.h>
#endif

#endif  //END GCTA_CPU_H