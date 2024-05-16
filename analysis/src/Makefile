# Makefile for building a shared library

# Compiler
# CXXFLAGS := $(shell root-config --cflags) #-ggdb 

CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,

LIB_SUFFIX=.so
LD = g++

# Fastjet
# FASTJET3_CONFIG=/lustre/isaac/scratch/tmengel/jet-background-subtraction/src/install/fastjet/bin
# CXX_FASTJET3 = $(shell $(FASTJET3_CONFIG)/fastjet-config --cxxflags --libs)

# ROOT
CXX_ROOT = $(shell root-config  --cflags --libs)

# Pythia8
# CXX_PYTHIA8 = $(shell pythia8-config --cxxflags --libs)

CXX_COMMON:= $(CXX_COMMON) $(GZIP_LIB)
# CXX_COMMON+= $(CXX_ROOT) $(CXX_PYTHIA8) $(CXX_FASTJET3)
# LIB_TOYMODEL = /lustre/isaac/scratch/tmengel/jet-vn-cumulant/external/RooUnfold/build
# CXX_COMMON += -I$(LIB_TOYMODEL) -L$(LIB_TOYMODEL) -Wl,-rpath,$(LIB_TOYMODEL) -lRooUnfold
CXX_COMMON+= $(CXX_ROOT)
# Source files

SRCS := CorrCalculator.cc \
		FlowCalculator.cc \
		MultSubCalibrator.cc \
		MultSubtractor.cc  \
		RefFlowCalculator.cc \
		ReweightJets.cc \
		SubSampler.cc \
		Unfolder.cc 

# Object files
OBJS := $(SRCS:.cc=.o)

# Shared library name
LIBNAME := libAna.so

# Target for building the shared library
$(LIBNAME): $(OBJS)
	$(CXX) -shared -o $@ $(OBJS) $(CXX_COMMON)

# Target for compiling individual source files
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $< -w $(CXX_COMMON)

# Target for cleaning up the build artifacts
clean:
	rm -f $(OBJS) $(LIBNAME)

# Phony target to avoid conflicts with files of the same name
.PHONY: clean
