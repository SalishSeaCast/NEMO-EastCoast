.. _Install_netcdfc:

*********************
Installing NETCDF-C
*********************

This is a sample script for installing NETCDF4-C using intel compilers. ::

 #Configuration options
 export CC=mpicc
 export CXX=mpicxx
 export FC=mpif90
 export F77=mpif90
 export F90=mpif90
 export LD=mpif90

 export CFLAGS='-O2 -fPIC '
 export CXXFLAGS="-O2 -fPIC "

 export F90FLAGS="-O2 -fPIC "
 export FCFLAGS="-O2 -fPIC "
 export FFLAGS="-O2 -fPIC "

 export LDFLAGS="-O2 -fPIC -shared-intel "
 # FLAGS FOR F90  TEST-EXAMPLES 
 export FCFLAGS_f90="-O2 -fPIC "

 #Installation Directories
 export HDF5DIR=home/fateme/Opt/opt-install/HDF5/HDF5_ZLIB_MPICH
 export NETCDFDIR=/home/fateme/Opt/opt-install/NETCDF4/NETCDF4_hdf5_mpich

 export CPPFLAGS="-I$HDF5DIR/include"
 export LDFLAGS=" -shared-intel -L$HDF5DIR/lib"

 export LD_LIBRARY_PATH=$HDF5DIR/lib:$LD_LIBRARY_PATH

 # Compile, test and install
 
 ./configure \
 --prefix=$NETCDFDIR \
 --disable-dap --disable-dap-remote-tests --enable-netcdf4 \
 --enable-shared --enable-parallel-tests  \
    2>&1 | tee fateme-configure.log

 make 2>&1 | tee fateme-make.log
 make check 2>&1 | tee fateme-check.log
 make install 2>&1 | tee fateme-install.log





