#this script was written by Dr. YAN Zeyin at 2019
#contact: zeyin.yan@outlook.com，Southern University of Science and Technology
#!/bin/bash
if [[ $# -le 0 ]]
then
	echo -e "\033[32mUsage: $0 xxx.gjf layer(H, M, L)\033[0m"
elif [[ $# -eq 2 ]]
then
	a=`sed -n '/[a-zA-Z0-9@#$%^&*]/!=' $1 | sed -n '2p'`
	b=`sed -n '/[a-zA-Z0-9@#$%^&*]/!=' $1 | sed -n '3p'`
	begin=$(($a+2))
	end=$(($b-1))
	line=$(($b-$begin))
	sed -n $begin','$end'p' $1 > coord.temp

	if [[ $2 == 'H' ]]
	then
		sed -e '/ L/d' coord.temp | sed -e '/ M/d' > char.temp
	elif [[ $2 == 'M' ]]
	then
		sed -e '/ L/d' coord.temp  > char.temp
	elif [[ $2 == 'L' ]]
	then
		mv coord.temp char.temp
	fi

	pos=`cat char.temp | awk '{print $1}' | awk -F '-' '{print $3}' | awk '{sum+=$1} END {print sum}'`
	neg=`cat char.temp | awk '{print $1}' | awk -F '--' '{print $2}' | awk '{sum+=$1} END {print sum}'`

	rm -f char.temp coord.temp
	charge=$( echo "$pos - $neg" | bc )
#	echo $pos
#	echo $neg
	echo 'The charge of layer: '$2 ' is ' $charge
fi



