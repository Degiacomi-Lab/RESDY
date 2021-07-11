#!/bin/bash
 	
# ...and collect results
FILES="*.pka"
for f in $FILES
do 
    n=$(echo $f | cut -d "." -f1)
    echo $n
    cat $f | grep "  LYS" | grep " 0 " | sed "s/LYS/$f LYS/" >> ../result.dat	
done


