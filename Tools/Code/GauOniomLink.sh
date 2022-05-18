#this script was written by Dr. YAN Zeyin at 2019
#contact: zeyin.yan@outlook.com，Southern University of Science and Technology
#!/bin/bash
if [[ $# -le 0 ]]
then
	echo -e "\033[32mUsage: $0 xxx.gjf layer(M, L)\033[0m"
elif [[ $# -eq 2 ]]
then
	a=`sed -n '/[a-zA-Z0-9@#$%^&*]/!=' $1 | sed -n '2p'`
	b=`sed -n '/[a-zA-Z0-9@#$%^&*]/!=' $1 | sed -n '3p'`
	begin=$(($a+2))
	end=$(($b-1))
	line=$(($b-$begin))

	if [[ $2 == 'M' ]]
	then
		line_num=`grep -n ' M H' $1 | awk -F ':' '{print $1}'`
	elif [[ $2 == 'L' ]]
	then
		line_num=`grep -n ' L H' $1 | awk -F ':' '{print $1}'`
	fi

	cp $1 $1_temp
	echo "This link atoms's line number are:"
	echo $line_num
	for x in $line_num
	{
		sed -i "$x,$x s/ -1 /  0 /g" $1_temp
	}
	echo 'difference:'
	diff $1 $1_temp

	sed -n $begin','$end'p' $1 | grep -n ' 0 ' | awk -F ':' '{print $1}' > relaxed0.temp
	num_line=`sed -n '$=' relaxed0.temp`
	echo $num_line > relaxed.temp
	cat relaxed0.temp >> relaxed.temp
fi



