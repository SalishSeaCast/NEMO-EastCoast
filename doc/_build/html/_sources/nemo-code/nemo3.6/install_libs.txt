******************************
Installing Required libraries
******************************

Required libraries
-------------------

To install NEMO v3.6, the following libraries are required:

* MPI
* HDF5
* NETCDF4
* XIOS

To compile NEMO3.6, XIOS IO server should be installed. XIOS needs NETCDF4 and if you want to use the "one_file" mode which means having one overall output instead of outputs for each processor, you will need the hdf/netcdf libraries properly compiled to allow parallel IO. In this document, sample scripts are provided for installing required libraries using Intel compilers.

Installing HDF5
---------------

Download the latest version of hdf5 from website:
http://www.hdfgroup.org/HDF5/release/obtain5.html#obtain

I downloaded version 1.8.13. Here is a sample script for installing this library: ::
 
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


Installing NETCDF4
------------------

To compile XIOS library we need NECTD4 and not NETCDF3. To Install NETCDF4, first download and install the developer version of netcdf-c from github:

https://github.com/Unidata/netcdf-c

As this version does not have a configuration file you first must make it by: ::

   autoreconf -i -f

The next step is to install netcdf-fortran. Download the latest stable version and install the library. A sample script is: ::
 



NOTE: Do not download the latest stable version of netcdf-c (4.3.2). If you do, you will encounter errors while compiling with enable-parallel option. This is due to a bug which has been fixed in the developers version.

Installing XIOS
---------------

Obtain the latest revision of XIOS: ::

   svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk XIOS

Follow instructions given here to install: http://forge.ipsl.jussieu.fr/ioserver/wiki/documentation


