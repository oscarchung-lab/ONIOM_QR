#/bin/bash
if [[ $# -eq 1 ]]
then
    gjff=$1
    name=${gjff%.*}
    a=`sed -n '/[a-zA-Z0-9@#$%^&*]/!=' $gjff | sed -n '2p'`
    b=`sed -n '/[a-zA-Z0-9@#$%^&*]/!=' $gjff | sed -n '3p'`
    begin=$(($a+2))
    end=$(($b-1))

    h2olineb=`grep -n '^[ ]*O-OW-' $gjff | head -1 | awk -F ":" '{print $1}'`
    h2olinee=`grep -n '^[ ]*H-HW-' $gjff | tail -1 | awk -F ":" '{print $1}'`

    h2ob=$(($h2olineb-$begin+1))
    h2oe=$(($h2olinee-$begin+1))

    c=`grep -n "^[ ]*$h2ob" $gjff | head -1 | awk -F ":" '{print $1}'`
    line=$(($c-1))


    head -n $line $gjff > $name-order.gjf 
    for ((i=$h2ob;i<$h2oe;i=i+3));
    do 
        numo=$i 
        numh1=$(($i+1))
        numh2=$(($i+2))
        echo '' $numo $numh1 '1.0' $numh2 '1.0' >> $name-order.gjf
        echo '' $numh1 >> $name-order.gjf
        echo '' $numh2 >> $name-order.gjf
    done
    echo >> $name-order.gjf
    echo >> $name-order.gjf
    #echo $h2olineb $h2olinee
    #echo $h2ob $h2oe
    #echo $c
else
    echo -e "\033[32mUsage: $0 xxx.gjf \033[0m"
fi