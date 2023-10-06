#!/bin/bash

#for ang in 1.0 1.25 1.5 1.538 1.58 1.75 2.0 3.0 4.0 5.0
for ang in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 10.0 20.0 30.0
do
   cd Ang$ang
   for i in `seq 0 1 100`
   do
   val=$(tail -1 force.${i}.txt)
   echo $val >> force.txt
   done
   cd ..
done
