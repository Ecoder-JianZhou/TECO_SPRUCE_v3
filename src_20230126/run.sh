export LD_LIBRARY_PATH=/home/jz964/miniconda3/lib/:$LD_LIBRARY_PATH 
gfortran -g -fbacktrace -Wall -fcheck=all dataType.f90 updateAndSummary.f90 writeOutputs2nc.f90 soil.f90 vegetation.f90 transfer.f90 driver.f90 mcmc.f90 spinup.f90  main.f90 -o run_teco -I/home/jz964/miniconda3/include -L/home/jz964/miniconda3/lib -lnetcdff -lnetcdf
# gfortran -g dataType.f90 driver.f90  mcmc.f90 spinup.f90 vegetation.f90 main.f90 -o run_teco
rm *.mod
./run_teco
rm run_teco
