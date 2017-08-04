# GCTA
Genome-wide Complex Trait Analysis (GCTA)

License: GPL v3.0

(C) 2010-2017, The University of Queensland

Please report bugs to: Jian Yang <jian.yang@uq.edu.au>

## **Develop page**

## Requirements
1. Only x86\_64 (Intel, or AMD64) based operation systems are supported. 
2. Intel MKL 2017 or above
3. EIGEN >= 3.3.3 
4. CMake >= 3.1
### Linux
1. Kernel version >= 2.6.8
2. GCC version >= 4.9 with C++ 11 support, Or Intel compiler version >= 15.0
### Mac
1. Mac version > 10.7
2. Mac platform SDK (Usually ship with Xcode)
3. GCC version >= 5.4 with C++ 11 support, Or Intel compiler version >= 15.0
### Windows
1. Windows 7 and above
2. MSbuild toolkit 2015 version >= 14.0 (Or Visual Studio 2015 and above)
3. Intel compiler version >= 15.0
4. zlib 


## Before compile
### Linux & Mac
1. clone the repository or download the code
2. revise environment variable  MKLROOT and EIGEN3\_INCLUDE\_DIR to your own location of MKL and EIGEN3

> Note: you can either export MKLROOT="LOCATION\_MKL" before compile, or put these variables permanently to shell profiles (.bashrc, .profile)

### Windows
1. clone the repository or download the code
2. revise the Include and Library path in gcta\_win64 folder

## Compile
### Linux & Mac
```
mkdir build
cd build
cmake ..
make
```
After it finishes, the binary file *gcta64* is in same folder

### Windows
MSbuild command line or Visual Studio build
After it finishes, the binary file is in the x64/Release folder

## Develop guide
There are many ways to contribute to GCTA, for minor changes, just change the code in the master branch. We'd better test that portion before commit back into the repository. 

For the libraries that use the partial loading manner, header files shall be put into *include* folder, and source code shall be put into *src* folder. CMake build system will recogonize them automatically in this way. Revision on the code of GCTA version 1.30 shall be in the *main* folder.

Don't change the directory structures, it will break the CMake configuration. If it is necessary, we'd better confirm everyone has pushed back the revisions, and we should change the CMake configuration accordingly. 

For huge changes, we can choose one of these methods: 
### Create new branch in the main repository
```
git branch branch_name
git checkout branch_name
```

> Notice: you'd better merge the changes from the master branch periodically. Otherwize, your branch is always built on old codes.

After it finishes, you'd better ask admin to merge it into the master branch if you are not confidence about the branch operations.

### Fork the repository by clicking the fork button
*Work in this way if you are new to git, as work directly on main repository may cause disaster*
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

2. The admin merge the pull request and merge the branch periodically to publish new version of GCTA. 


> Notice: pay attention to the changes in the main master repository,
> We should always fecth the upstream and merge the changes first before making the **New pull request**,
> or the admin may face lots of conflictions.
> Always make **New pull request** when the function can really work after full testing.
