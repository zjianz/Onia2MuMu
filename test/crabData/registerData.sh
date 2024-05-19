#!/bin/bash
# A script to generate crab config file automaticly

# project config
dataList='Run2023dataList.txt'
template='crab3_template.py'
fileName='crab3'
cat $dataList | while read rows
do
	# tag1,2,3 are tags from CMS data naming standard
	# tag1: physics target, tag2: data set, tag3: data format
	tag1=$(echo $rows | awk 'BEGIN{FS="/"} {print $2}')
	tag2=$(echo $rows | awk 'BEGIN{FS="/"} {print $3}')
	tag3=$(echo $rows | awk 'BEGIN{FS="/"} {print $4}')
	# define tag format with varibale tag, it will show up at TaskTag in the config .py files and also in their filenames 
	tag1_=$(echo $tag1 | grep -o [0-9]) 
	tag2_=$(echo $tag2 | awk 'BEGIN{FS="-"} {print $1$3}')
	tag=${tag1_}_${tag2_}_$tag3

	sed -e 's:DataSet:'"$rows:" $template > ${fileName}_${tag}.py
	sed -i -e 's:TaskTag:'"${fileName}_${tag}:" ${fileName}_${tag}.py
	#mv ${fileName}_${tag}.py result/
done
