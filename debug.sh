#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
cd $BASEDIR
which cmake
cmake --version 2>&1
tmp=$?
if [ ${tmp} -ne 0 ]
then
    echo "ERROR: You should install cmake(2.8 or later) first."
    echo "  Please goto https://cmake.org to download and install it."
    exit
fi

rm -rf build
mkdir build
cd build

cmake .. -DCMAKE_BUILD_TYPE=Debug

tmp=$?
echo "cmake compile return:" ${tmp}
if [ ${tmp} -ne 0 ]
then
    echo "cmake compile <>return:" ${tmp}
    exit -1
fi

make
tmp=$?
echo "make return:" ${tmp}
if [ ${tmp} -ne 0 ]
then
    echo "make <>return:" ${tmp}
    exit -1
fi

cd ../check/bin
valgrind --leak-check=full ./netrec "MODE=3 CUST=1000 SUPP=100 BUDGET=50 TYPE=1 PROB1=0.1 PROB2=0.9 SAMP=10 QOS=10 TIME=300"
