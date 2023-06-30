#those steps tested on Ubuntu 18. and all linux distribution should works.

cd ~
mkdir gcta_build
cd gcta_build
rootdir=`pwd`

cd $rootdir
wget https://yanglab.westlake.edu.cn/data/gcta_dep.tar
tar -xf gcta_dep.tar

cd $rootdir/gcta_dep
#if you have problems to do this please just copy the mkl file to here $rootdir/mkl_pkg
./l_onemkl_p_2022.0.1.117_offline.sh -a -c --install-dir $rootdir/mkl_pkg

cd $rootdir/gcta_dep
tar -zxf boost_1_75_0.tar.gz
cd boost_1_75_0
./bootstrap.sh --prefix=$rootdir/boost_pkg
./b2
./b2 install

cd $rootdir/gcta_dep
tar -zxf eigen-3.3.7.tar.gz
cd eigen-3.3.7
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$rootdir/eigen_pkg ..
make
make install

cd $rootdir/gcta_dep
tar -zxf Spectra_v1.0.0.tar.gz
cd spectra-1.0.0/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$rootdir/spectra_pkg ..
make
make install

cd $rootdir/gcta_dep
tar -zxf gsl-2.7.tar.gz
cd gsl-2.7/
./configure --prefix=$rootdir/gsl_pkg
make
make install
export CPATH="$rootdir/gsl_pkg/include":$CPATH
export LIBRARY_PATH="$rootdir/gsl_pkg/lib":$LIBRARY_PATH
export LD_LIBRARY_PATH="$rootdir/gsl_pkg/lib":$LD_LIBRARY_PATH

cd $rootdir/gcta_dep
tar -zxf zlib-1.2.11.tar.gz
cd zlib-1.2.11/
./configure  --prefix=$rootdir/zlib_pkg
make
make install
export CPATH="$rootdir/zlib_pkg/include":$CPATH
export LIBRARY_PATH="$rootdir/zlib_pkg/lib":$LIBRARY_PATH
export LD_LIBRARY_PATH="$rootdir/zlib_pkg/lib":$LD_LIBRARY_PATH

cd $rootdir/gcta_dep
tar -zxf zstd_v1.5.0.tar.gz
cd zstd-1.5.0/build/cmake/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX="$rootdir/zstd_pkg" ..
make
make install
export CPATH="$rootdir/zstd_pkg/include":$CPATH
export LIBRARY_PATH="$rootdir/zstd_pkg/lib":$LIBRARY_PATH
export LD_LIBRARY_PATH="$rootdir/zstd_pkg/lib":$LD_LIBRARY_PATH

cd $rootdir/gcta_dep
tar -zxf tcl8.6.11-src.tar.gz
cd tcl8.6.11/unix/
./configure --prefix=$rootdir/tcl_pkg
make
make install
export PATH=$rootdir/tcl_pkg/bin:$PATH
export CPATH="$rootdir/tcl_pkg/include":$CPATH
export LIBRARY_PATH="$rootdir/tcl_pkg/lib":$LIBRARY_PATH

cd $rootdir/gcta_dep
tar -zxf sqlite.tar.gz
cd sqlite
./configure --prefix=$rootdir/sqlite_pkg
make
make install
export CPATH="$rootdir/sqlite_pkg/include":$CPATH
export LIBRARY_PATH="$rootdir/sqlite_pkg/lib":$LIBRARY_PATH
export LD_LIBRARY_PATH="$rootdir/sqlite_pkg/lib":$LD_LIBRARY_PATH

cd $rootdir
git clone https://github.com/jianyangqt/gcta.git
cd gcta
git submodule update --init
export EIGEN3_INCLUDE_DIR=$rootdir/eigen_pkg/include/eigen3
export SPECTRA_LIB=$rootdir/spectra_pkg/include
export BOOST_LIB=$rootdir/boost_pkg/include
export MKLROOT=$rootdir/mkl_pkg/mkl/latest
mkdir build
cd build
cmake ..
make
