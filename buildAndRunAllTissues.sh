#!/bin/bash
rm -rf OutputFiles/
rm -rf build/
rm -rf Tissue*/
mkdir build
cmake -Bbuild -Hsrc
cd build
make
for i in `seq 1 76`;
do
mkdir ../OutputFiles
#cp ../Recurrence/tissuePar_$i.dat ../InputFiles/tissuePar.dat
#cp ../Recurrence/treatment_$i.dat ../InputFiles/treatment.dat
#cp ../HistSpec/tissueDim$i.dat ../InputFiles/tissueDim.dat
#cp ../HistSpec/inTum$i.dat ../InputFiles/inTum.dat
#cp ../HistSpec/inVes$i.dat ../InputFiles/inVes.dat
#cp ../HistSpec/inPO2$i.dat ../InputFiles/inPO2.dat
cp ../Recurrence/tissueParADCT2w_$i.dat ../InputFiles/tissuePar.dat
./RTsim

cp ../InputFiles/tissueDimRec.dat ../Recurrence/tissueDimRec_$i.dat
cp ../InputFiles/inTumRec.dat ../Recurrence/inTumRec_$i.dat
cp ../InputFiles/inVesRec.dat ../Recurrence/inVesRec_$i.dat
mv ../OutputFiles ../Tissue$i
done

