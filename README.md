# GCTA
Genome-wide Complex Trait Analysis (GCTA)


(C) 2010-Present, Jian Yang <jian.yang.qt@gmail.com>

2016-Present, Zhili Zheng <zhilizheng@outlook.com> (New codes from version 1.90)

Please report bugs to: Jian Yang <jian.yang.qt@gmail.com>

License: GPLv3 (default license); some codes are released under LGPL detailed in files.

## **Develop page**

## Requirements
1. Only x86\_64 (Intel, or AMD64) based operation systems are supported. 
2. Intel MKL 2017 or above
3. Eigen == 3.3.7 (bugs in new version eigen)
4. CMake >= 3.1
5. BOOST >= 1.4
6. zlib >= 1.2.11
7. sqlite3 >= 3.31.1 
8. zstd >= 1.4.4
9. Spectra >= 0.8.1

### Linux
1. Kernel version >= 2.6.28, or the MKL library can't work. 
2. GCC version >= 6.1 with C++ 11 support.

### Mac & Windows
We don't provide the instructions to compile on these platforms due to the complexity and huge version difference. Please download the binary from our website.

## Before compile
Environment variable (Linux only) 
1. point MKLROOT and EIGEN3\_INCLUDE\_DIR to your own location of MKL and EIGEN3
2. point BOOST\_LIB to your own location of BOOST 
3. point SPECTRA\_LIB to the path of Spectra

```
# this will update the submodule
# If you download the code archieve:
cd SOURCE_CODE_FOLDER_DECOMPRESSED
git clone https://github.com/zhilizheng/plink-ng submods/plink-ng

# If you clone from github
git clone THIS_REPO_ADDRESS
git submodule update --init

# make the folder
mkdir build
cd build
cmake ..
make
```
After it finishes, the binary file *gcta64* is located in same folder

## Develop guide
### Fork the repository by clicking the fork button
*Work in this way if you are new to git, as working directly on main repository may cause disaster*
1. Revise the code 
    * Minor bug fix
    1. Clone your own forked repository.
    2. Revise the code in the master branch, and push back into your own repository.
    3. Make a **New pull request** from your own repository, write the detailed log of changes, why revise it, how to use it
    
    ---

    * A new function
    1. Clone your own forked repository.
    2. New branch, write the code. Push back into your own repository.
    3. Make a **New pull request** from your own repository, write the document how to use the new function

2. I will merge the pull request and merge the branch periodically to publish new version of GCTA. 


> Note: pay attention to the changes in the main master repository,
> We should always fecth the upstream and merge the changes first before making the **New pull request**,
> or the admin may meet lots of conflicts.
> Always make **New pull request** when the function can really work after full testing.
