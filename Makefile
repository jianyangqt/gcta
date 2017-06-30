# -----------------------------------------------------------------
#   Makefile for GCTA 
#   
#   Supported platforms (x86_64 only)
#       Linux                       UNIX
#       Windows                     WIN
#       MacOS                       MAC
# ---------------------------------------------------------------------

# change the MKL ROOT before compile, we encourage to use 2017 version and above
# MKLROOT = $(ib)/local/intel_MKL2017_U3/compilers_and_libraries/linux/mkl
MKLROOT = /opt/intel/compilers_and_libraries/mac/mkl

# change the EIGEN path before compile, we encourage to use 3.3.3 version and above
#EIGEN = $(ib)/local/Eigen/eigen-3.3.4
EIGEN = /usr/local/Cellar/eigen/3.3.4/include/eigen3

# Use sinlge precision to store matrix
#SINGLE_PRECISION = 1 
ifdef SINGLE_PRECISION
 CXXFLAGS += -DSINGLE_PRECISION=1
endif

OUTPUT = ./release/gcta64

############################################
###  Linux configuration ##################

CXXFLAGS = -w -s -O3 -m64 -fopenmp -DNDEBUG -msse2 -std=c++11 -DMKL_LP64 -I. -I$(MKLROOT)/include -I$(EIGEN) 
LDFLAGS = -static -lz -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl

############################################
### Mac and Windows configuration ##########
### Windows not supported by this makefile # 

# Windows specified 
ifeq ($(OS),Windows_NT)
    # zlib header
    # CXXFLAGS += -I zlib_header
    # LDFLAGS = mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib zlib.lib
    $(error Windows is not supported by make, you can turn to gcta_win64 to build)
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        # cann't compile with clang, we must install g++ manually
        CXX = g++-7
	# static linking have some problems
        LDFLAGS = -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L/usr/local/opt/llvm/lib -liomp5 -lpthread -lm -ldl -lz 
    endif
endif

############################################
### Configuration of all platform ##########


HDR += CommFunc.h \
	   cdflib.h \
	   dcdflib.h \
	   eigen_func.h \
           gcta.h \
	   ipmpar.h \
           StatFunc.h \
           StrFunc.h \
           zfstream.h
SRC = bivar_reml.cpp \
           CommFunc.cpp \
	   eigen_func.cpp \
           data.cpp \
	   dcdflib.cpp \
           edata.cpp \
           ejma.cpp \
           est_hsq.cpp \
           grm.cpp \
           gwas_simu.cpp \
           ld.cpp \
           joint_meta.cpp \
           mlm_assoc.cpp \
	   mkl.cpp \
           option.cpp \
           popu_genet.cpp \
           raw_geno.cpp \
           sbat.cpp \
           StatFunc.cpp \
           StrFunc.cpp \
           reml_within_family.cpp \
           zfstream.cpp
	   
OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT) 

$(OUTPUT) :
	mkdir -p release
	$(CXX) -o $(OUTPUT) $(OBJ) $(LDFLAGS) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean:
	rm -f *.o *~
