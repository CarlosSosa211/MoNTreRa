path=~/Carlos/Results/Recurrence/simp/TTum330_alphaG1120_vd025
for k in `seq 0 4`
do
./buildAndRunAllTissues.sh
mkdir $path/rep$k
cp -rf Tissue* $path/rep$k/
for i in `seq 1 76`
do
for file in $path/rep$k/Tissue$i/*
do 
mv "$file" "${file/.res/_$i.res}"
done
cp $path/rep$k/Tissue$i/* $path/rep$k/
done
rm -rf $path/rep$k/Tissue*
done


