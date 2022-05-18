#!/bin/bash
#Desc: Run PDB2PQR
#mode1： pdb
#mode2： pdb + pH
#mode3： pdb + mol2
#mode4： pdb + pH + mol2
source /share/env/python-3.6.8.env

if [[ $# -eq 1 ]]
then 
    mode=1 #pdb
elif [[ $# -eq 2 ]]
then 
    #if grep '^[[:digit:]]*$' <<< `sed 's/.//g' $2` >> /dev/null
    var2=${2:0-4}
    #echo $var2
    #if [[ "$var2"x = 'mol2'x ]]
    if [[ "$var2" == 'mol2' ]]
    then 
        mode=3  #pdb + mol2
    else
        mode=2 #pdb + pH
    fi

elif [[ $# -eq 3 ]]
then 
    mode=4  #pdb + pH + mol2
else
    echo -e "\033[32mUsage: $0 PDBfile \033[0m"
    echo -e "\033[32mUsage: $0 PDBfile pH\033[0m"
    echo -e "\033[32mUsage: $0 PDBfile mol2\033[0m"
    echo -e "\033[32mUsage: $0 PDBfile pH mol2\033[0m"
    echo -e "\033[32mUsage:PQB2PQR: default using AMBER at pH of crystallization\033[0m"
    #echo -e "\033[32mUsage:PQB2PQR: default using AMBER at pH of crystallization\033[0m"
    exit 101
fi

echo $mode

fpdb=$1
fpqr=`echo $1 | sed 's/.pdb/.pqr/g'`
if [[ mode -eq 1 || mode -eq 3 ]]
then
    ph=`grep 'CRYSTALLIZATION CONDITIONS:' $fpdb | awk -F "PH" '{print $2}' | awk '{print $1}' | sed 's/,//g'`
    
    if [ -z $ph ]
    then
        ph='7.0'
    fi
else
    ph=$2
fi
echo 'PH= '$ph

#run
if [[ mode -eq 1 || mode -eq 2 ]]
then
    pdb2pqr30 --ff=AMBER --ffout=AMBER --titration-state-method=propka --with-ph=$ph $fpdb $fpqr 
else
    fmol=${!#}
    if [[ -f "$fmol" ]]
    then
        pdb2pqr30 --ff=AMBER --ffout=AMBER --ligand=$fmol --titration-state-method=propka --with-ph=$ph $fpdb $fpqr
    else
        echo -e "\033[32m $fmol does not exist, Run with not ligand file\033[0m"
        pdb2pqr30 --ff=AMBER --ffout=AMBER --titration-state-method=propka --with-ph=$ph $fpdb $fpqr
    fi 
fi
#--ligand=LIGAND_FILE   Assign charges to the ligand specified in a MOL2 file

