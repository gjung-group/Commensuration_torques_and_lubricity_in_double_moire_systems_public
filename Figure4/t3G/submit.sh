#!/bin/sh

#for ang in 0.1 #0.5 1.0 5.0 10.0 20.0 30.0
#for ang in 0.5 1.0 5.0 10.0 20.0 30.0
#for ang in 1.0 
#for ang in 0.0 1.0 1.25 1.5 1.538 1.58 1.75 2.0 3.0 4.0 5.0
#for ang in 1.875 2.5 3.5
#for ang in 0.25 0.5 0.75
for ang in 0.1
#for ang in 1.5385 #0 0.5 1 1.5385 2 2.3  #  0.5 1 #1.56 #4.0 5.0
do
   mkdir -p Ang$ang
   cd Ang$ang
   cp -p ../in.quasistatic .
   cp -p ../C* .
   cp -p ../B* .
   cp -p ../r.sh .
   cp -p ../Gen* .
   cp -p ../mol2xyz.py .
   sed -i "s/XXX/${ang}/g" GenMoire2Flake4.f90
   gfortran GenMoire2Flake4.f90
   ./a.out > cellInfo
   {
   tail -1 BLBL.mol
   } && {
   #grep -A3 "Final lattices" cellInfo
   grep -A4 "Final N atoms" cellInfo
           #python fracToCart.py BLBL.frac
   #python replaceFrac.py
           mv BLBL.mol BLBLNew.mol
   val0=$(grep 'Total' cellInfo |tr -d ' '| cut -d ":" -f 2)
   sed -i "s/NNN/${val0}/g" BLBLNew.mol
   python mol2xyz.py
 
 
   AAA=$(grep 'FixA:' cellInfo | tr -d ' '| cut -d ',' -f 2) 
   BBB=$(grep 'FixA:' cellInfo | tr -d ' '| cut -d ',' -f 3) 
   CCC=$(grep 'FixA:' cellInfo | tr -d ' '| cut -d ',' -f 4) 
   DDD=$(grep 'FixA:' cellInfo | tr -d ' '| cut -d ',' -f 5) 
 
   FFF=$(grep 'FixD' cellInfo | tr -d ' '| cut -d ',' -f 2) # FixD : (a1+a2)/30  ~ 0.14 Angstrom 
   GGG=$(grep 'FixD' cellInfo | tr -d ' '| cut -d ',' -f 3) 
 
#   FFF=$(grep 'FixC' cellInfo | tr -d ' '| cut -d ',' -f 2) # FixC : (a1+a2)/210 ~ 0.02 Angstrom 
#   GGG=$(grep 'FixC' cellInfo | tr -d ' '| cut -d ',' -f 3) 
  
   NNN=100

   sed -i "s/AAA/${AAA}/g"  in.quasistatic
   sed -i "s/BBB/${BBB}/g"  in.quasistatic
   sed -i "s/CCC/${CCC}/g"  in.quasistatic
   sed -i "s/DDD/${DDD}/g"  in.quasistatic
   sed -i "s/EEE/${val0}/g" in.quasistatic
   sed -i "s/FFF/${FFF}/g"  in.quasistatic
   sed -i "s/GGG/${GGG}/g"  in.quasistatic
   sed -i "s/NNN/${NNN}/g"  in.quasistatic
   }
   sbatch r.sh
   cd ..
done
