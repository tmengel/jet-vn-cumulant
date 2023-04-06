# Compilation flags (change these for your system)
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=


FASTJET3_CONFIG=/lustre/isaac/scratch/tmengel/jet-background-subtraction/src/install/fastjet/bin/fastjet-config
FASTJET3_BIN=
FASTJET3_INCLUDE=-I/lustre/isaac/scratch/tmengel/jet-background-subtraction/src/install/fastjet/include
FASTJET3_LIB=-L/lustre/isaac/scratch/tmengel/jet-background-subtraction/src/install/fastjet/lib -Wl,-rpath,/lustre/isaac/scratch/tmengel/jet-background-subtraction/src/install/fastjet/lib -lfastjet

ROOT_CONFIG=/nfs/home/tmengel/.conda/envs/myhepenv/bin/root-config
ROOT_BIN=/nfs/home/tmengel/.conda/envs/myhepenv/bin/
ROOT_INCLUDE=-I/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/include
ROOT_LIB=-L/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/lib -Wl,-rpath,/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/lib -lCore

PYTHIA8_CONFIG=/nfs/home/tmengel/.conda/envs/myhepenv/bin/pythia8-config
PYTHIA8_BIN=/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/bin/
PYTHIA8_INCLUDE=-I/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/include
PYTHIA8_LIB=-L/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/lib -Wl,-rpath,/lustre/isaac/scratch/tmengel/.conda/envs/myhepenv/lib
