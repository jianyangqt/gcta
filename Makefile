
# -----------------------------------------------------------------
#   Makefile for GCTA 
#   
#   Supported platforms
#       Unix / Linux                UNIX
#       Windows                     WIN
#       MacOS                       MAC
# ---------------------------------------------------------------------

# Set this variable to either UNIX, MAC or WIN
SYS = UNIX
OUTPUT = ./gcta64_testing

MKLROOT = /opt/intel/mkl

# Use sinlge precision to store matrix
SINGLE_PRECISION = 

# Put C++ compiler here; Windows has it's own specific version
CXX_UNIX = g++
CXX_WIN = C:\CodeBlocks\MinGW\bin\mingw32-g++.exe
CXX_MAC = g++

# Any other compiler flags here ( -Wall, -g, etc)
CXXFLAGS = -w -O3 -m64 -fopenmp -I ~/Lib/eigen -DEIGEN_NO_DEBUG -msse2 -std=c++0x -I.

ifdef SINGLE_PRECISION
 CXXFLAGS += -DSINGLE_PRECISION=1
endif

# Some system specific flags

ifeq ($(SYS),WIN)
 CXXFLAGS += -DWIN -static -I ~/Lib/zlib
 LIB += ~/Lib/zlib/zlib.lib 
 CXX = $(CXX_WIN)
endif

ifeq ($(SYS),UNIX)
 CXXFLAGS += -DUNIX -static -m64 -I$(MKLROOT)/include
 LIB += -static -lz -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
 CXX = $(CXX_UNIX)
endif

ifeq ($(SYS),MAC)
 CXXFLAGS += -DUNIX -Dfopen64=fopen 
#-isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5
 LIB += -lz
 CXX = $(CXX_MAC)
endif

#ifeq ($(SYS),SOLARIS)
# CXX = $(CXX_UNIX)
#endif

HDR += CommFunc.h \
	   cdflib.h \
	   dcdflib.h \
           gcta.h \
	   ipmpar.h \
           SimplexSolver.h \
           StatFunc.h \
           StrFunc.h \
           zfstream.h
SRC = bivar_reml.cpp \
           CommFunc.cpp \
           data.cpp \
	   dcdflib.cpp \
           edata.cpp \
           ejma.cpp \
           est_hsq.cpp \
           grm.cpp \
           grm_wt.cpp \
           gwas_simu.cpp \
           ld.cpp \
           joint_meta.cpp \
           mlm_assoc.cpp \
	   mkl.cpp \
           option.cpp \
           popu_genet.cpp \
           raw_geno.cpp \
           sbat.cpp \
           SimplexSolver.cpp \
           StatFunc.cpp \
           StrFunc.cpp \
           reml_within_family.cpp \
           zfstream.cpp
	   
OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT) 

$(OUTPUT) :
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJ) $(LIB) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean:
	rm -f *.o *~
