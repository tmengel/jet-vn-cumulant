################################################################################
## Makefile for jet-cumulants/simulation/src
################################################################################
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,

LIB_SUFFIX=.so
LD = g++

# Fastjet
FASTJET3_CONFIG=/lustre/isaac/scratch/tmengel/jet-background-subtraction/src/install/fastjet/bin
CXX_FASTJET3 = $(shell $(FASTJET3_CONFIG)/fastjet-config --cxxflags --libs)

# ROOT
CXX_ROOT = $(shell root-config  --cflags --libs)

# Pythia8
CXX_PYTHIA8 = $(shell pythia8-config --cxxflags --libs)

# ToyModel
LIB_TOYMODEL=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/jet-vn-cumulant/simulation/src 
CXX_TOYMODEL = -I$(LIB_TOYMODEL) -L$(LIB_TOYMODEL) -Wl,-rpath,$(LIB_TOYMODEL) -lToyModel

LIB_ANA=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/jet-vn-cumulant/analysis/src
CXX_ANA = -I$(LIB_ANA) -L$(LIB_ANA) -Wl,-rpath,$(LIB_ANA) -lAna

# CXX_COMMON:= $(CXX_COMMON) $(GZIP_LIB) $(CXX_ROOT)
# CXX_SIM:= $(CXX_COMMON) $(GZIP_LIB) $(CXX_PYTHIA8) $(CXX_FASTJET3) $(CXX_TOYMODEL)
# CXX_ANA:= $(CXX_COMMON) $(CXX_ANA)
CXX_COMMON:= $(CXX_COMMON) $(CXX_ROOT) $(CXX_PYTHIA8) $(CXX_FASTJET3) $(CXX_TOYMODEL) $(CXX_ANA)
################################################################################

# Targets
BINS = $(patsubst %.cpp, %, $(wildcard *.cpp))
	
default : $(BINS)

$(BINS): % : %.cpp
	@echo -n "Building $@ ... "
	$(CXX) $<  -o $@ -w  $(CXX_COMMON)
	@echo "Done"

clean: 
	rm -f *.o	$(BINS)

jobclean:
	rm *.o* *.e*
