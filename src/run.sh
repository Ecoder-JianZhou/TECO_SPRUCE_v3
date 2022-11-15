# gfortran -g dataType.f90 driver.f90  mcmc.f90 spinup.f90 driver.f90 vegetation.f90 soil.f90 transfer.f90 main.f90 -o run_teco
gfortran -g dataType.f90 driver.f90  mcmc.f90 spinup.f90 vegetation.f90 main.f90 -o run_teco
./run_teco
