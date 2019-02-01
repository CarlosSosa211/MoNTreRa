#!/bin/bash
rm -rf OutputFiles/
rm -rf build/
mkdir build
cmake -Bbuild -Hsrc
cd build
make
for i in `seq 1 21`;
do
mkdir ../OutputFiles
cp ../HistSpec/tissueDim$i.dat ../InputFiles/tissueDim.dat
cp ../HistSpec/inTum$i.dat ../InputFiles/inTum.dat
cp ../HistSpec/inVes$i.dat ../InputFiles/inVes.dat
./RTsim
mv ../OutputFiles ../Tissue$i
done

