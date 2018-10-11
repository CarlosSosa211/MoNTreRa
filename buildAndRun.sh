#!/bin/bash
rm -rf build/
mkdir build
rm -rf OutputFiles
mkdir OutputFiles
cmake -Bbuild -Hsrc
cd build
make
./RTsim

