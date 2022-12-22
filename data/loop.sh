#!/bin/bash

FILES=$(ls | grep data-rough)

for VAR in $FILES
do
   VARNAME="${VAR%.*}"
   ./check.py $VAR
   #./rough-distr.py $VAR > "$VARNAME-rough.txt"
done
