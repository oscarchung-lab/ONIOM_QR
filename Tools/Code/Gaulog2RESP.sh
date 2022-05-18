#/bin/bash

if [[ $# -eq 0 ]]
then
    echo -e "\033[32mUsage: $0 Gau.log\033[0m"
    exit 101
else
    fname=`echo $1 | sed 's/.log//g'`
    NC=`grep ' Charge =' $1 | awk '{print $3}'`
    M=`grep ' Charge =' $1 | awk '{print $6}'`
    antechamber -i $1 -fi gout -o $fname.mol2 -fo mol2 -c resp -nc $NC -m $M -rn MOL
    parmchk -i $fname.mol2 -f mol2 -o $fname.frcmod
fi

