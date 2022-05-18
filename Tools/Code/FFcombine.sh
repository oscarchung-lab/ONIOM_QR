#!/bin/bash

if [[ $# -eq 2 ]]
then
    #split files to diff parts
    grep '^[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]*' $1 > atom.tmp1
    grep '^[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]*' $1 > bond.tmp1
    grep '^[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]*' $1 > angle.tmp1
    grep '^[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9\-][0-9]*\.[0-9]*' $1 > torsion.tmp1
    grep '^[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9\-][0-9]*\.[0-9]*' $1 > improper.tmp1
    grep '^  [A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]*' $1 > vdw.tmp1


    grep '^[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]*' $2 > atom.tmp2
    grep '^[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]*' $2 > bond.tmp2
    grep '^[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]*' $2 > angle.tmp2
    grep '^[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9\-][0-9]*\.[0-9]*' $2 > torsion.tmp2
    grep '^[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ]-[A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9\-][0-9]*\.[0-9]*' $2 > improper.tmp2
    grep '^  [A-Z0-9][A-Za-z0-9\* ] [ ]*[0-9][0-9]*\.[0-9]* [ ]*[0-9][0-9]*\.[0-9]*' $2 > vdw.tmp2

    :>CobineFF.dat 
    echo $1 '+' $2 'combined parameter' >> CobineFF.dat 
    #atoms
    echo 'MASS' >> CobineFF.dat 
    while read line
    do
        stratom=${line:0:8}
        yn=`grep "^$stratom" atom.tmp2`
        if [  -z "$yn" ]
        then 
            printf "$line \n" >> CobineFF.dat
        else 
            printf "  $line \n"
        fi
    done  < atom.tmp1
    cat atom.tmp2 >> CobineFF.dat
    echo >> CobineFF.dat

    #bonds
    echo 'BOND' >> CobineFF.dat 
    while read line
    do
        strbond=${line:0:5}
        yn=`grep "^$strbond" bond.tmp2`
        if [  -z "$yn" ]
        then 
            printf "$line \n" >> CobineFF.dat
        else 
            printf "  $line \n"
        fi
    done  < bond.tmp1
    cat bond.tmp2 >> CobineFF.dat
    echo >> CobineFF.dat

    #angles
    echo 'ANGL' >> CobineFF.dat 
    while read line
    do
        strangle=${line:0:8}
        yn=`grep "^$strangle" angle.tmp2`
        if [  -z "$yn" ]
        then 
            printf "$line \n" >> CobineFF.dat
        else 
            printf "  $line \n"
        fi
    done  < angle.tmp1
    cat angle.tmp2 >> CobineFF.dat
    echo >> CobineFF.dat

    #torsion
    echo 'torsion' >> CobineFF.dat 
    while read line
    do
        strtorsion=${line:0:11}
        yn=`grep "^$strtorsion" torsion.tmp2`
        if [  -z "$yn" ]
        then 
            printf "$line \n" >> CobineFF.dat
        else 
            printf "  $line \n"
        fi
    done  < torsion.tmp1
    cat torsion.tmp2 >> CobineFF.dat
    echo >> CobineFF.dat

    #improper
    echo 'improper' >> CobineFF.dat 
    while read line
    do
        strimproper=${line:0:11}
        yn=`grep "^$strimproper" improper.tmp2`
        if [  -z "$yn" ]
        then 
            printf "$line \n" >> CobineFF.dat
        else 
            printf "  $line \n"
        fi
    done  < improper.tmp1
    cat improper.tmp2 >> CobineFF.dat
    echo >> CobineFF.dat

    #vDw
    echo 'vDw' >> CobineFF.dat 
    while read line
    do
        strvdw=${line:0:4}
        yn=`grep "^$strvdw" vdw.tmp2`
        if [  -z "$yn" ]
        then 
            printf "  $line \n" >> CobineFF.dat
        else 
            printf "  $line \n"
        fi
    done  < vdw.tmp1
    cat vdw.tmp2 >> CobineFF.dat
    echo >> CobineFF.dat

    rm -f *.tmp1 *.tmp2 
else
    echo -e "\033[32mUsage: $0 parm1 parm2 \033[0m"
    echo -e "\033[33mif conflicts, the parm2 is favorable \033[0m"
    exit 101
fi
