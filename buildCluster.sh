#!/bin/bash
sed -i "s/tissueDim.*\.dat/tissueDim$1.dat/g" src/main.cpp
sed -i "s/inTum.*\.dat/inTum$1.dat/g" src/main.cpp
sed -i "s/inVes.*\.dat/inVes$1.dat/g" src/main.cpp
sed -i "s/_[0-9][0-9]*\.res/_$1.res/g" src/morris.cpp
cmake -Bbuild -Hsrc
cd build
make
mv RTsim RTsim$1
