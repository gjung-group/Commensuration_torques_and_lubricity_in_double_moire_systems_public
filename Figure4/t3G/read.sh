#!/bin/bash

#for ang in 0.0 1.0 1.25 1.5 1.538 1.58 1.75 2.0 3.0 4.0
#for ang in 1.875 2.5 3.5
for ang in 0.0 0.1 0.25 0.5 0.75 1.5 1.58 2.25 1.125 2.75 1.0 1.25 1.5 1.538 1.58 1.75 2.0 3.0 4.0 1.875 2.5 3.5 2.25 1.125 2.75
#for ang in 0.0
#for ang in 1.538
do
   cd Ang$ang
   echo $ang
   rm force.txt
   #for i in `seq 0 1 45`
   for i in `seq 0 1 15`
   do
   val=$(tail -1 force.${i}.txt)
   echo $val >> force.txt
   done
   cd ..
done
