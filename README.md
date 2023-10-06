# Commensuration torques and lubricity in double moire systems

This folder contains scripts that can be used to reproduce the results from the paper **Commensuration torques and lubricity in double moire systems** available [here](https://arxiv.org/pdf/2301.04105.pdf).

The files are organized by figure for which they are used. We provide here further details on how to use them and/or refer to the comments provided in the files.

## Figure 1
- The `LAMMPS` calculations can be performed using the files provided in `mirrorSymmetryBrokent3G` subfolder. One neds to relax in two consecutive steps (`relaxAll` is the second step after the first step breaks the mirror symmetry by fixing the bottom layer) to obtain the final structure. Relaxation conditions are given in the two `lammps_cg.in` input files. The different `python` scripts can be used to extract the output from the first run and give it as input for the next run.
- Panels (a) and (b) were prepared using VESTA and/or Keynote.
- Panel (c) and panel (d) can be obtained by using the python script `getInterlayerDistance.py` where interpolation methods are used to get the z-coordinate of each layer and the vertical distance between two layers.

## Figure 2
- We performed a lot of `LAMMPS` calculations for this figure. The indices that control the commensurate simulations for these calculations can be found in the Appendix from our paper. These indices can be defined directly in the `GenMoire2.f90` file (replacing the values in the current `n = (/ 58, 58, 58, 60, 59, 59 /)` line) or one can comment this line and use the `FindMoire` and `FindMoire2` routines to find the commensurate cells for any lattice mismach and angle combination as defined by the hard-coded values of `degtol`, `lattol`, `phideg`, `aG`, `aBN0`, as well as the counterparts for the second moire.
- To control the rotation center, we fix three corner atoms of the periodic cell to keep their x-y coordinates in the right stacking. This is controlled by the `fix` flags in the `lammps_cg.in` file.
- For the sliding maps, we vary the values of the `basis2` in the `GenMoire.f90` for the second to be able to map all possible stackings.
- For the energy cuts from panel (c), we provide a series of scripts and a jupyter notebook in the `energyContributions` subfolder. This process involves several caerful pre-processing and post-processing steps, so we refer to the paper the understand how each contribution is defined.

## Figure 3
We refer to the paper and reference therein for the definition of how one atom is assigned to a specific high-symmetry stacking family. This definition requires the extraction of the displacement vectors from the relaxed positions. This can be achieved using the example given in the `stackingAreasNotebook.ipynb` jupyter notebook which requires the initial and final position `generate.xyz` and `generateInit.xyz` files extracted from a `LAMMPS` calculation using `getOutput.py` (see Figure 1 script for reference).

## Figure 4
- To run the `LAMMPS` quasi-static friction calculations, one should run the `submit.sh` file in either the `t2G` ir `t3G` subfolders modifying the script to include the angles one is interested in. These values will automatically be parseed to obtain the commensurate cells using the `GenMoire2Flake4.f90` file. 
- Once the calculations are ready, one can use `read.sh` to obtain the force.txt files for each angle. 
- These can then be further processed and visualized using the `finalNotebookForPaperFigures.ipynb` jupyter notebook.
- Figure preparation was done in the keynote file provided as `finalFigurePreparation.key`
