# load modules
module load legacy
module load intel
module load hdf5

# export environmental variables
export OMP_NUM_THREADS=2

export I2STM_LIBS="-L/opt/intel/Compiler/11.1/075/lib/intel64 -lmkl_n1ntel_lp64 -lmkl_n1ntel_thread -lmklhore -lpthread -lifcore -mcmodel=medium -openmp -lhdf5"

# run
#bsub -I -W 00:10 -n 4 -R 'rusage[mem=6144]' -J rsfhk ./i2evps

bsub -W 6:00 -n 4 -R 'rusage[mem=10000,scratch=2400]' -J VEP-t1 ./i2evps

#bsub -W 24:00 -n 4 -R 'rusage[mem=10000,scratch=2400]' -J VEP-t1 ./i2evps
#bsub -W 24:00 -n 4 -R 'rusage[mem=10000,scratch=2400]' -J VEP-t1 -w "ended(VEP-t1)" ./i2evps
#bsub -W 24:00 -n 4 -R 'rusage[mem=10000,scratch=2400]' -J VEP-t1 -w "ended(VEP-t1)" ./i2evps
#bsub -W 24:00 -n 4 -R 'rusage[mem=10000,scratch=2400]' -J VEP-t1 -w "ended(VEP-t1)" ./i2evps

echo ">>> submission of jobs..."
