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

cp ../Recurrence/tissueDimRec3DReal_$i.dat ../InputFiles/tissueDim.dat
cp ../Recurrence/inTumRec3DReal_$i.dat ../InputFiles/inTumDim.dat
cp ../Recurrence/inVesRec3DReal_$i.dat ../InputFiles/inVesDim.dat
cp ../Recurrence/treatment_$i.dat ../InputFiles/treatment.dat
./RTsim
mv ../OutputFiles ../Tissue$i

#cp ../Recurrence/tissueDimRec3DReal_$i.dat ../InputFiles/tissueDim.dat
#cp ../Recurrence/inTumRec3DReal_$i.dat ../InputFiles/inTum.dat
#./RTsim
#cp ../InputFiles/inVes.dat ../Recurrence/inVesRec3DReal_$i.dat


#cp ../Recurrence/tissueParSpheADCT2w_$i.dat ../InputFiles/tissuePar.dat
#cp ../Recurrence/treatment_$i.dat ../InputFiles/treatment.dat
#cp ../HistSpec/tissueDim$i.dat ../InputFiles/tissueDim.dat
#cp ../HistSpec/inTum$i.dat ../InputFiles/inTum.dat
#cp ../HistSpec/inVes$i.dat ../InputFiles/inVes.dat
#cp ../HistSpec/inPO2$i.dat ../InputFiles/inPO2.dat
#./RTsim

#cp ../InputFiles/tissueDim.dat ../Recurrence/tissueDimRec_$i.dat
#cp ../InputFiles/inTum.dat ../Recurrence/inTumRec_$i.dat
#cp ../InputFiles/inVes.dat ../Recurrence/inVesRec_$i.dat
#mv ../OutputFiles ../Tissue$i
done

