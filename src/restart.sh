# load modules
module load legacy
module load intel
module load hdf5

# export environmental variables
export OMP_NUM_THREADS=2

export I2STM_LIBS="-L/opt/intel/Compiler/11.1/075/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lifcore -mcmodel=medium -openmp -lhdf5"

# delete previous output files
# rm -rf Maxslipvel.txt 
rm -rf *lsf*

echo "=================================== "
echo ">>> Compile in2"
icc -std=c99 in2evps.c -o in2evps -ltiff ${I2STM_LIBS}
echo ">>> ...OK!"
echo ">>> Compile i2"
icc -std=c99 i2evps.c -o i2evps -ltiff ${I2STM_LIBS}
echo ">>> ...OK!"
echo "=================================== "
