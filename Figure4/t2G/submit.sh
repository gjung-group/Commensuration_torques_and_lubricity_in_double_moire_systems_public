#!/bin/sh

#for ang in 0.1 #0.5 1.0 5.0 10.0 20.0 30.0
#for ang in 0.0 0.5 1.0 5.0 10.0 20.0 30.0
#for ang in 0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9
for ang in 1.5 2.0 2.5 3.0 3.5 4.0 4.5
do
   mkdir -p Ang$ang
   cd Ang$ang
   cp ../in.* .
   cp ../C* .
   cp ../r.sh .
   cp ../Gen* .
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
   }
   sbatch r.sh
   cd ..
done
