# VELO2CYCLEs
Visco-ElastO 2d CYCLe of Earthquakes


Computational Earthquake Physics
ETH Zurich, 2022
Dal Zilio, L., Lapusta, N., Avouac, J. P., Gerya, T. (2022). 
Subduction earthquake sequences in a non-linear visco-elasto-plastic megathrust. 
Geophysical Journal International, 229(2), 1098-1121.
DOI: https://doi.org/10.1093/gji/ggab521

![examples](https://github.com/lucadalzilio/VELO2CYCLEs/blob/main/banner/velo2cycles_banner.png)

=====================================================

This software uses the following libraries:
```
module load intel 
module load hdf5  
```
Additionally this software uses the Pardiso solver that uses MKL.
Be aware that you need a valid intel licence to use those.

=====================================================

The c version of this simulation works on the ETH EULER cluster only.
To run it on a local device or another cluster, the .sh scripts would need to be changed accordingly.
  
some important intel flags:
```
  –lmkl_intel_thread 
  -lmkl_core -lpthread 
  -lifcore 
  -mcmodel=medium 
  -qopenmp 
  -lhdf5
```

=====================================================

To run the code on ETH cluster EULER use the

and the .sh script provided in the corresponding folder of this project. 
(These scripts are written for the EULER cluster only!)

```
source compile_euler.sh
source run.sh
```
=====================================================
