#!/bin/bash
OPTIND=$1
# set default values to 0
if [ -z $OPTIND ]; then
    OPTIND=0
fi

SRCDIR=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/jet-vn-cumulant/analysis/src
MACRODIR=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/macros/ana 
if [ $OPTIND -eq 1 ]; then
    echo "Making Ana"
    cd $SRCDIR
    make clean
    make -j4

    cd $MACRODIR
    make clean
    make -j4
else
    echo "Making Ana"
    cd $SRCDIR
    make -j4

    cd $MACRODIR
    make -j4
fi