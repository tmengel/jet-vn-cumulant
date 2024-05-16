################################################################################
## Makefile for jet-cumulants/simulation/src
################################################################################
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,

LIB_SUFFIX=.so
LD = g++

# ROOT
CXX_ROOT = $(shell root-config  --cflags --libs)

# Pythia8

CXX_COMMON:= $(CXX_COMMON) $(GZIP_LIB)
# LIB_ROOUNFOLD=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/external/RooUnfold/build
# CXX_COMMON += -I$(LIB_ROOUNFOLD) -L$(LIB_ROOUNFOLD) -Wl,-rpath,$(LIB_ROOUNFOLD) -lRooUnfold
CXX_COMMON+= $(CXX_ROOT)
LIB_ANA=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/ana/src
CXX_COMMON += -I$(LIB_ANA) -L$(LIB_ANA) -Wl,-rpath,$(LIB_ANA) -lAna



################################################################################

# Targets
BINS = $(patsubst %.cpp, %, $(wildcard *.cpp))
	
default : $(BINS)

$(BINS): % : %.cpp
	@echo -n "Building $@ ... "
	$(CXX) $(CXXFLAGS) $< -o  $@ -w  $(CXX_COMMON)
	@echo "Done"

clean: 
	rm -f *.o	$(BINS)

jobclean:
	rm *.o* *.e*
