export LD_LIBRARY_PATH=/home/jz964/miniconda3/lib/:$LD_LIBRARY_PATH 
gfortran newTest_WriteTime.f90 -o newtest -I/home/jz964/miniconda3/include -L/home/jz964/miniconda3/lib -lnetcdff -lnetcdf
./newtest
