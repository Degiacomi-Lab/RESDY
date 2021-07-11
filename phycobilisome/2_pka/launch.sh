#!/bin/bash

# launch calculations...
FILES="../selections/*pdb"
for f in $FILES
do
  echo "Processing $f file..."

  # swap methylated asparagine for asparagine (extra methyl already removed in VMD)
  sed $f 's/MEN/ASN/g'

  fname=$(echo $f | cut -d "/" -f3)
  n=$(echo $fname | cut -d "." -f1)
  logfile=$n".log"
  pkafile=$n".pka"

  propka3 $f > $logfile
  cat $pkafile | grep "  LYS" | grep " 0 " | sed "s/LYS/$n LYS/" | awk 'BEGIN { OFS = "\t"} {print $1, $2, $3, $5}' >> ../result.dat

done


# checks
ls *log | xargs grep " 0 " | grep -v "LYS" | cut -d "." -f1 | uniq > problem_minor.dat
ls *log | xargs grep " 0 " | grep "LYS" | cut -d "." -f1 | uniq > problem_major.dat
