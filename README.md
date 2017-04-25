# GCTA

Develop page

## Requirements
1. Linux x64
2. GCC version >= 4.9 
3. MKL 2017 (Old version have issue when inverse large matrix)
4. EIGEN >= 3.3.3

## Before compile
1. clone the repository
2. revise the Makefile, specify the path of MKLROOT and EIGEN

## Compile
```
make
```
After it finishes, the binary file are in release folder

## Develop guide
1. Fork the repository into your own account by clicking the fork button
2. Revise the code
    * Minor bug fix
    1. Clone your own forked repository.
    2. Revise the code in the master branch, and push back into your own repository.
    3. Make a **New pull request** from your own repository, write the detailed log of changes, why revise it, how to use it

    * A new function
    1. Clone your own forked repository.
    2. New branch, write the code. Push back into your own repository.
    3. Make a **New pull request** from your own repository, write the document how to use the new function

3. The admin merge the pull request and merge the branch periodically to publish new version of GCTA. 

> Notice: pay attention to the changes in the main master repository,
we should always fecth the upstream and merge the changes first before making the **New pull request**,
or the admin will face lots of conflictions.

> Always make **New pull request** when the function can really work after full testing.
