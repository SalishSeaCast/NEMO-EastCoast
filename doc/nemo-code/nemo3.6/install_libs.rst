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

I downloaded version 1.8.13. A sample script for installing this library can be found in :ref:`Install_hdf5`

Installing NETCDF4
------------------

To compile XIOS library we need NECTD4 and not NETCDF3. To Install NETCDF4, first download and install the developer version of netcdf-c from github:

https://github.com/Unidata/netcdf-c

As this version does not have a configuration file you first must make it by: ::

   autoreconf -i -f

A sample script for installing netcdf-c is here: :ref:`Install_netcdfc`

The next step is to install netcdf-fortran. Download the latest stable version and install the library. A sample script can be found in :ref:`Install_netcdff`

NOTE: Do not download the latest stable version of netcdf-c (4.3.2). If you do, you will encounter errors while compiling with enable-parallel option. This is due to a bug which has been fixed in the developers version.

Installing XIOS
---------------

Obtain the latest revision of XIOS: ::

   svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk XIOS

Follow instructions given here to install: http://forge.ipsl.jussieu.fr/ioserver/wiki/documentation


