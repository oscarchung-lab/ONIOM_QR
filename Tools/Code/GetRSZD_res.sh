#/bin/bash

if [[ $# -eq 0 ]]
then
    echo -e "\033[32mUsage: $0 PDB including all key residues\033[0m"
    exit 101
else
    awk '($1=="ATOM" || $1=="HETATM") {print $0}' $1 | cut -b18-20 > resname.tmp
    awk '($1=="ATOM" || $1=="HETATM") {print $0}' $1 | cut -b23-26 > resnum.tmp
    paste resname.tmp resnum.tmp | uniq | sed 's/ //g' > resinfo.tmp
    #paste resname.tmp resnum.tmp | uniq | sed 's/ //g' | sed 's/\t//g'> resinfo.tmp
    numres=`cat resinfo.tmp | wc -l`

    echo 'methods' > title.tmp 
    echo $numres 'residues'
    for logf in *.txt
    do 
        echo $logf
        echo "${logf%.*}" >> title.tmp
        :>scores.tmp
        for ((i=1; i<=$numres; i++))
        do 
            #echo $i
            resn=`sed -n $i'p' resinfo.tmp | awk '{print $1}'`
            resnum=`sed -n $i'p' resinfo.tmp | awk '{print $2}'`
            #echo $resn $resnum
            awk '($1=='\"$resn\"' && $3=='\"$resnum\"') {print $(NF-5)}' $logf >> scores.tmp

            
            
        done
        paste resinfo.tmp scores.tmp | sed 's/ //g' > resinfo.tmp1
        mv resinfo.tmp1 resinfo.tmp
    done 
fi

tr "\n" " " < title.tmp > RSZDScores.dat 
echo >> RSZDScores.dat 
cat resinfo.tmp >> RSZDScores.dat 

rm -rf *.tmp 
