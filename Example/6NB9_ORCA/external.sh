#!/bin/bash
#orca as QM program
inpf='6nb9.inp'
#Update external QM program input file (geometry info from coord.xyz)
jobn=`echo $inpf | sed 's/.inp//g'`
atomn=`head -n 1 coord.xyz | tr -cd "[0-9]"`
a=`grep -rin xyz $inpf | awk -F ":" '{print $1}'`
b=$(($atomn+2))

head -n $a $inpf > temp.inp
sed -n '3,'$b'p' coord.xyz >> temp.inp
echo '*' >> temp.inp
mv temp.inp $inpf

cp $PBS_O_WORKDIR/*.inp $orca_dir/ 2>/dev/null
cp $PBS_O_WORKDIR/*.gbw $orca_dir/ 2>/dev/null
cp $PBS_O_WORKDIR/*.xyz $orca_dir/ 2>/dev/null

cd $orca_dir
ls -l > $PBS_O_WORKDIR/orca.log
echo $orca_dir >> $PBS_O_WORKDIR/orca.log
#call external QM program



$ORCA_BIN/orca $jobn.inp >> $PBS_O_WORKDIR/orca.log

cp $orca_dir/*.gbw $PBS_O_WORKDIR/ 2>/dev/null
cp $orca_dir/*.xyz $PBS_O_WORKDIR/ 2>/dev/null
cp $orca_dir/*.engrad $PBS_O_WORKDIR/ 2>/dev/null

cd $PBS_O_WORKDIR
#output the energy and gradient results in $jobn.engrad to the QM_grad.dat file

sed -n '8p' $jobn.engrad > QM_grad.dat
sed -n '12,/^#/p' $jobn.engrad | sed '$d' | awk 'ORS=NR%3?"\t":"\n"' >> QM_grad.dat
