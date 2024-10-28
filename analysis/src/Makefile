# Makefile for building a shared library
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,

LIB_SUFFIX=.so
LD = g++

# ROOT
CXX_ROOT = $(shell root-config  --cflags --libs)

CXX_COMMON:= $(CXX_COMMON) $(GZIP_LIB)
CXX_COMMON+= $(CXX_ROOT)

SRCS := CorrCalculator.cc \
		DiffCorrFromHistos.cc \
		FlowCalculator.cc \
		MultSubCalibrator.cc \
		MultSubtractor.cc  \
		RefFlowCalculator.cc \
		ReweightJets.cc \
		SimpleDiffFlow.cc \
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
