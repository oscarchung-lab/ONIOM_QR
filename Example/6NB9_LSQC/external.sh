#!/bin/bash
#lsqc as QM program
#define the parameter of the gjf file
#atom=224
title=10
#xyz='2ona.xyz'
gjf="6nb9.gjf"
atom=`head -n 1 coord.xyz | tr -cd "[0-9]" `
#calculate number
a=$(($title+1))
b=$(($atom+$a))
c=$(($atom+3))
dir=${gjf%.*}
#Update external QM program input file (geometry info from coord.xyz)
if [ ! -d "$dir" ]; then
    head -n $title $gjf > temp.gjf
else
    cp $dir/$dir/$dir.cha .
    cp $dir/$dir/$dir.frg .
    cp $dir/$dir.gjf .
    rm -rf $dir
    head -n $title $gjf | sed s/'\scharge=[a-zA-Z]*'/' charge=read'/g | sed s/'\sfrag=[a-zA-Z]*'/' frag=read'/g > temp.gjf
fi
sed -n '3,'$c'p' coord.xyz >> temp.gjf
sed -n $b',/^$/p' $gjf >> temp.gjf
mv temp.gjf $gjf

#call external QM program
lsqc $gjf > $dir.out

#output the energy and gradient results in 2ona.engrad to the QM_grad.dat file
cp $dir/$dir/$dir.force .
cp $dir/$dir/$dir.out .
head -n 1 $dir.out | awk '{print $4}' > QM_grad.dat
force_grad $dir.force $atom
cat $dir.force >> QM_grad.dat
