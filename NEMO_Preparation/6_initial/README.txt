NAME: gen_ini.m

 AUTHOR: J.-P. Paquin 

 DATE: Feb 2014 

 REVISIONS: 

 DESCRIPTION: Generate Initial conditions for Regional NEMO ("DESTINATION"). 
              Developped to take data from GLORYS version 2v3 Data  ("SOURCE")
              from ${datastorage}/DATA/GLORYS2V3
              - One script for all variables T,S,U,V,SSH

 HYPOTHESES:  


 INPUTS:  INPUT dataset to extract the initial conditions (IC)
            OPTION1: coordinates of SOURCE grid
            OPTION2: files containing latitudes and longitudes for
                     U, V, T and F on SOURCE grid 
            DATA COVERAGE : 1993-2011 (monthly or annual)

            DESTINATION 'NEMO regional' grid: coordinates, bathymetry, mesh_mask 


 OUTPUTS: Interpolated initial conditions for U and V on DESTINATION grid 


 CALLED PGM & SCRIPTS: f_open_netcdf
                       f_readnetcdf
                       f_writenetcdf
                       f_readmeshmask

                    PROCESSING T, S, and SSH
                       f_create_mask (for floodnan_opa3)
                       f_interp_scalar (floodnan_opa3; interp1q)

                    PROCESSING U and V
                       f_interpuv (->f_rot_rep->f_opa_angle)


 NOTES: REQUIRES M_MAP PACKAGE (version 1.4) for Matlab

