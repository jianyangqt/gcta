# GCTA
GCTA (Genome-wide Complex Trait Analysis) is a software package, which was initially developed to estimate the proportion of phenotypic variance explained by all genome-wide SNPs for a complex trait but has been extensively extended for many other analyses of data from genome-wide association studies (GWASs). Please see the software website through the link below for more information.

Software website: https://yanglab.westlake.edu.cn/software/gcta/
License: GPLv3 (some parts of the code are released under LGPL as detailed in the files).


## Credits  
Jian Yang developed the original version (before v1.90) of the software (with supports from Peter Visscher, Mike Goddard and Hong Lee) and currently maintains the software.

Zhili Zheng programmed the fastGWA, fastGWA-GLMM and fastGWA-BB modules, rewrote the I/O and GRM modules, improved the GREML and bivariate GREML modules, extended the PCA module, and improved the SBLUP module.  

Zhihong Zhu programmed the mtCOJO and GSMR modules and improved the COJO module.  

Longda Jiang and Hailing Fang developed the ACAT-V module.  

Jian Zeng rewrote the GCTA-HEreg module.  

Andrew Bakshi contributed to the GCTA-fastBAT module.

Angli Xue improved the GSMR module.

Robert Maier improved the GCTA-SBLUP module.

Contributions to the development of the methods implemented in GCTA (e.g., GREML methods, COJO, mtCOJO, MLMA-LOCO, fastBAT, fastGWA and fastGWA-GLMM) can be found in the corresponding publications (https://yanglab.westlake.edu.cn/software/gcta/index.html#Overview).


## Questions and Help Requests
If you have any bug reports or questions please send an email to Jian Yang at <jian.yang@westlake.edu.cn>.


## Compilation

#### Requirements
1. Currently only x86\_64-based operating systems are supported.
2. Intel MKL 2017 or above
3. Eigen == 3.3.7 (there are bugs in the new version of Eigen)
4. CMake >= 3.1
5. BOOST >= 1.4
6. zlib >= 1.2.11
7. sqlite3 >= 3.31.1
8. zstd >= 1.4.4
9. Spectra >= 0.8.1
10. gsl (GNU scientific library)

#### Linux
1. Kernel version >= 2.6.28 (otherwise the Intel MKL library doesn't work).
2. GCC version >= 6.1 with C++ 11 support.

#### Mac OS & Windows
Here we do not provide instructions to compile GCTA in Mac OS or Windows because of the enormous complexity of the compilation process and differences among OS versions. We suggest you download the compiled executable files directly from the GCTA website.

#### Before compilation (Linux)
1. Export MKLROOT variable into shell environment  
`export MKLROOT=where_your_MKL_located`    
2. Export EIGENE3_INCLUDE_DIR into your shell environment  
`export EIGENE3_INCLUDE_DIR=where_your_eigen3_located`  
3. Export BOOST_LIB environment variable  
`export BOOST_LIB=where_your_boost_include_directory_located`  
4. Export SPECTRA_LIB environment variable
`export SPECTRA_LIB=where_your_spectra_lib_located`
5. Other compilation Requirements should be placed in your compiler's header files and library searching paths.  
6. Clone GCTA source code from github 
`git clone https://github.com/jianyangqt/gcta GCTA_PATH`
7. Clone plink_ng submods
```
cd GCTA_PATH/submods
git clone https://github.com/zhilizheng/plink-ng
```

#### To compile
```
cd GCTA_PATH
mkdir build
cd build
cmake ..
make
```
