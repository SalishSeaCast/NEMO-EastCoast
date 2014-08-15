************************************
Compiling NEMOv3.6
************************************

This is a document on compiling NEMO 3.6. It also includes which versions of required libraries are appropriate for this task. Example scripts for compiling the libraries are also provided. These scripts are based on Mercator-ocean's scripts for compiling NEMOv3.5.

Getting the code
=================

Download the latest version from the trunk repository: ::

    svn --username yourusername co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM

Your username is the same as the one you use for http://www.nemo-ocean.eu/ website


Required libraries
==================

To install NEMO v3.6, the following libraries are required:

* HDF5
* NETCDF4
* XIOS

To compile NEMO3.6, XIOS IO server should be installed. XIOS need NETCDF4 and if you want to use the "one_file" mode which means having one overall output instead of outputs for each processor, you will need the hdf/netcdf libraries properly compiled to allow parallel IO. 

Installing Required libraries
=============================

Installing HDF5
---------------

Download the latest version of hdf5 from website:
http://www.hdfgroup.org/HDF5/release/obtain5.html#obtain

I downloaded version 1.8.13. You can find the script for installing HDF5 here:

Installing NETCDF4
------------------

To compile XIOS library we need NECTD4 and not NETCDF3. To Install NETCDF4, first download and install the developer version of netcdf-c from github:

https://github.com/Unidata/netcdf-c

As this version does not have a configuration file you first must make it by: ::

   autoreconf -i -f

The script for installing netcdf-c library can be found here:

The next step is to install netcdf-fortran. Download the latest stable version and install the library. A sample script for this step can be found here:

NOTE: Do not download the latest stable version of netcdf-c (4.3.2). If you do, you will encounter errors while compiling with enable-parallel option. This is due to a bug which has been fixed in the developers version.

Installing XIOS
---------------

Obtain the latest revision of XIOS: ::

   svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk XIOS

Follow instructions given here to install: http://forge.ipsl.jussieu.fr/ioserver/wiki/documentation

Compiling NEMOv3.6
------------------
To compile NEMOv3.6 you first need to create an architecture file compatible with your machine which also indicates the path to NETCDF4 and XIOS libraries. You can find example architecture files in NEMOGCM/ARCH folder. 

After creating the architecture file, compile and create executable using an existing configuration. For example to use GYRE configuration and create a configuration called MY_GYRE: ::
  cd NEMOGCM/CONFIG
  ./makenemo –m your_architecture –r GYRE -n MY_GYRE

Please refer to http://www.nemo-ocean.eu/Using-NEMO/User-Guides/Basics/NEMO-Quick-Start-Guide for further details





