#!/bin/bash
#xtb as QM program
#Update external QM program input file (geometry info from coord.xyz)
#coord.xyz can be directly used

#call external QM program
#xtb has been added to the PATH
xtb coord.xyz --chrg 1 --uhf 0 --grad > xtb.out

#output the energy and gradient results in 2ona.engrad to the QM_grad.dat file
atom=`head -n 1 coord.xyz | tr -cd "[0-9]"`

#echo $atom

a=$(($atom+3))
b=$(($atom+$atom+2))

#echo $a $b # gradient file will be written continuously the next time
#rm -f gradient
#tail -n $b gradient | grep 'SCF energy' | awk '{print $7}' > QM_grad.dat
#tail -n $b gradient | sed -n $a','$c'p' >> QM_grad.dat
grep 'SCF energy' gradient | awk '{print $7}' > QM_grad.dat
sed -n $a','$b'p' gradient >> QM_grad.dat


rm -f charges energy xtbrestart gradient hessian xtbout mol.xyz tmpxx vibspectrum
rm -f hessian xtb_normalmodes g98_canmode.out g98.out wbo xtbhess.coord .tmpxtbmodef

