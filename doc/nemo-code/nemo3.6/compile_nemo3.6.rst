************************************
Compiling NEMOv3.6
************************************

Getting the code
=================

Download the latest version from the trunk repository: ::

    svn --username yourusername co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM

Your username is the same as the one you use for http://www.nemo-ocean.eu/ website


Compiling code
==================
To compile NEMOv3.6 you first need to create an architecture file compatible with your machine which also indicates the path to NETCDF4 and XIOS libraries. You can find example architecture files in NEMOGCM/ARCH folder. 

After creating the architecture file, compile and create executable using an existing configuration. 
For example to use GYRE configuration and create a configuration called MY_GYRE: ::
   
   cd NEMOGCM/CONFIG
   ./makenemo –m your_architecture –r GYRE -n MY_GYRE

Please refer to http://www.nemo-ocean.eu/Using-NEMO/User-Guides/Basics/NEMO-Quick-Start-Guide for further details





