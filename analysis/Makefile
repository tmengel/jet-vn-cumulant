# MakeFile
-include ../Makefile.inc

# CXX=g++
# CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
# CXX_SHARED=-shared
# CXX_SONAME=-Wl,-soname,
# LIB_SUFFIX=.so
LD = g++
LDFLAGS = 

CCFLAGS = $(shell root-config --cflags) #-ggdb 

# FASTJET3_CONFIG=/lustre/isaac/scratch/tmengel/jet-background-subtraction/src/install/fastjet/bin
# FASTJET3_INCLUDE=-I/lustre/isaac/scratch/tmengel/jet-background-subtraction/src/install/fastjet/include
# FASTJET3_LIB=-L/lustre/isaac/scratch/tmengel/jet-background-subtraction/src/install/fastjet/lib -Wl,-rpath,/lustre/isaac/scratch/tmengel/jet-background-subtraction/src/install/fastjet/lib -lfastjet

# ROOT_CONFIG=/nfs/home/tmengel/.conda/envs/myhepenv/bin/root-config
# ROOT_BIN=/nfs/home/tmengel/.conda/envs/myhepenv/bin/
# ROOT_INCLUDE=-I/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/include
# ROOT_LIB=-L/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/lib -Wl,-rpath,/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/lib -lCore

# PYTHIA8_CONFIG=/nfs/home/tmengel/.conda/envs/myhepenv/bin/pythia8-config
# PYTHIA8_BIN=/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/bin/
# PYTHIA8_INCLUDE=-I/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/include
# PYTHIA8_LIB=-L/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/lib -Wl,-rpath,/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/lib


BINS = $(patsubst %.C, %, $(wildcard *.C))

CXX_PYTHIA= $(PYTHIA8_INCLUDE) $(PYTHIA8_LIB) -lpythia8 -ldl

LIBS = $(shell root-config --libs) -lMLP -lMinuit -lTreePlayer -lTMVA -lTMVAGui -lXMLIO  -lMLP -lm
CXX_ROOT = $(ROOT_LIB) $(shell root-config --cflags --glibs)
default : $(BINS)

$(BINS): % : %.C 
	@echo -n "Building $@ ... "
	$(CXX)  $(CCFLAGS) $<  -o $@ -w  $(CXX_COMMON) $(CXX_PYTHIA) $(CXX_ROOT) $(FASTJET3_INCLUDE) $(FASTJET3_LIB) $(shell $(FASTJET3_CONFIG)/fastjet-config --cxxflags --libs)
	@echo "Done"


clean: 
	rm -f *.o	$(BINS)

jobclean:
	rm *.o* *.e*
