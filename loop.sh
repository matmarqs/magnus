#!/bin/bash

FILES=$(ls data | grep data-rough.txt)

make clean && make

for VAR in $FILES
do
   echo $VAR
   VARNAME="${VAR%.*}"
   ./main data/"$VAR" > output/Magnus4-"$VARNAME"_dotprod.txt
   #./check.py $VAR
   #./rough-distr.py $VAR > "$VARNAME-rough.txt"
done
