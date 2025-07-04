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

Wujuan Zhong and Judong Shen programmed the fastGWA-GE module. 

Contributions to the development of the methods implemented in GCTA (e.g., GREML methods, COJO, mtCOJO, MLMA-LOCO, fastBAT, fastGWA and fastGWA-GLMM) can be found in the corresponding publications (https://yanglab.westlake.edu.cn/software/gcta/index.html#Overview).


## Questions and Help Requests
If you have any bug reports or questions please send an email to Jian Yang at <jian.yang@westlake.edu.cn>.

## Code
```bash
# Clone gcta along with it's submodules
git clone --recursive https://github.com/yywei/gcta.git
```

## Dependencies
You can use a simple [script](https://github.com/yywei/gcta/blob/master/script/install_prerequisites.sh) to install them, which supports the following package managers: apt, vcpkg, and brew.
```bash
# See what packages will be install
./scripts/install_prerequisites.sh -l

# Override the package manager choice and install packages
# if one has not already been selected manually
./scripts/install_prerequisites.sh -m apt
```
You'll see the dependencies are
- cmake
- intel mkl (x64)
- openblas (Apple Silicon)
- boost
- eigen
- zlib
- zstd
- gsl
- sqlite3
## Building

#### Linux
```bash
# get code
cd ~/your_directory
git clone --recursive https://github.com/yywei/gcta.git
cd gcta

# Install dependencies
./scripts/install_prerequisites.sh

#initializing oneAPI environment
source /opt/intel/oneapi/setvars.sh

# Configure and build
cmake -B build
cmake --build build
```

#### MacOS
```bash
# get code
cd ~/your_directory
git clone --recursive https://github.com/yywei/gcta.git
cd gcta

# Install dependencies
./scripts/install_prerequisites.sh

# Configure and build
# Not support AppleClang
cmake -B build \
    -D $(brew --prefix llvm)/bin/clang \
    -D $(brew --prefix llvm)/bin/clang++
cmake --build build
```

#### Windwos
On Windows, you need to install [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) and [CMake](https://cmake.org/download/) separately. We recommended Git Bash for executing the bash snippets on this page. On Windwos, we use ***vcpkg*** package manager to install other dependencies, you should ensure that the directory containing ***vcpkg.exe*** has been added to the ***Path*** environment variable. Please note that MinGW and MSVC are not supported; only Clang/LLVM is. We recommend using Visual Studio 2022, and you can find instructions on how to install Clang/LLVM [here](https://learn.microsoft.com/en-us/cpp/build/clang-support-msbuild?view=msvc-170).
```bash
# get code
cd ~/your_directory
git clone --recursive https://github.com/yywei/gcta.git
cd gcta

# Install vcpkg
./script/vcpkg/bootstrap-vcpkg.sh

# Ensure that the directory containing vcpkg.exe has been added to the Path environment variable.

# Install dependencies
./scripts/install_prerequisites.sh

# Configure and build
cmake -B build \
    -T ClangCL \
    -G "Visual Studio 17 2022" \
    -D CMAKE_TOOLCHAIN_FILE=~/your_vcpkg_directory/scripts/buildsystems/vcpkg.cmake
cmake --build build --config Release
```
