#!/bin/bash

cmake --version

bld=r
BUILD_TYPE=Release

if [ "$1" == "d" ] ; then
    bld=d
    BUILD_TYPE=Debug
fi

mkdir -p build

cd build

cmake .. -DCMAKE_BUILD_TYPE=$BUILD_TYPE

cmake --build .