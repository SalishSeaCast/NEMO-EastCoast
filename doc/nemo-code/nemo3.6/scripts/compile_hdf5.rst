
Installing HDF5
================

This is a sample script for installing HDF5 using intel compilers. ::

 #!/bin/bash

 export HDF5_Make_Ignore=yes
 export CC=mpicc
 export CXX=mpicxx
 export FC=mpif90

 # Configure

 ./configure --prefix=/home/fateme/Opt/opt-install/HDF5/HDF5_ZLIB_MPICH \
 --enable-fortran  --enable-parallel --enable-hl --enable-shared  \
 2>&1 | tee fateme-configure_hdf5_zlib.log

 # Make and install

 make 2>&1 | tee fateme-make_hdf5_zlib.log

 make check  2>&1 | tee fateme-make_check_hdf5_zlib.log

 make install 2>&1 | tee fateme-install_hdf5_zlib.log

