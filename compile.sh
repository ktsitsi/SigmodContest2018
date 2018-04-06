#!/usr/bin/env bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
TYPE="Release"

cd $DIR
mkdir -p build/"${TYPE,,}"
cd build/"${TYPE,,}"
cmake -DCMAKE_BUILD_TYPE=${TYPE} -DFORCE_TESTS=OFF ../../src
make -j8