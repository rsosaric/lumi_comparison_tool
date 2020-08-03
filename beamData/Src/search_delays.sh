#!/bin/bash

awk -F "\"*,\"*" -v del=$1 'BEGIN{Tmp=0; Tmp2=0} NR>2{if (($2 != 0) && (Tmp2 != 0) && ($5 - Tmp > del)) print Tmp, $5, ($5 - Tmp)} NR>1{Tmp=$5; Tmp2=$2}' $2

