#!/bin/bash
#=======================================================
# Description :
#    Script for processing Herman Bathymetry Data 
#
# Authors: 
#    Ji Lei (Ji.Lei@dfo-mpo.gc.ca)                
#    Fatemeh Chegini (fatemeh.chegini@dal.ca)
#
# File Input:
#    NABD (Hermann) Data file. 
#    Each line of the file should have 14 columns
#    The header of the file is:
#    Year Month Day Hr Min Sec Lon  ...
#    Lat Obj_type Med_depth Agency  ...
#    Depth_type Vertical_ref Group_Id
#    Note: by default Med_depth is negative down
#
# Namelist Input:
#   *Selected Domain Area: 
#      minimum and maximum longitude and latitude
#   *Name of Input Herman Data file
#   *Outputdir : Directory to store your output files
#   *These inputs are read from 0_input.sh 
#   *0_input.sh is created when running 
#   Create_Bathymetry_fromNABDandEtopo1_Main.m script
# 
# Processing Steps:
#   1. Extract data that are in Selected Domain 
#   2. Eliminate bad and interpolated data
#     Eliminate lines in file that: 
#     do not have 14 columns
#     do not have proper reference level
#     have suspisous unit for depth 
#          i.e. depth is greater than 10,000 m       
#     Eliminate all interpolated data
#   3. Separate data that should be corrected with WebTide 
#     if reference level is mean sea level 
#     or depth > 200m data do not need to correction
#     if reference level is LLWLT 
#     and depth < 200m data should be corrected
#
#====================================================================

. ./0_Input.sh

#=========== Name of Output files which are created ==================

ExtractedDomain=$outputdir'/ExtractedDomainData.dat' 
GoodData=$outputdir'/GoodDatainDomain.dat'
BadData=$outputdir'/BadDatainDomain.dat'
MSLData=$outputdir'/MSLData.lld'
LLWLTData=$outputdir'/LLWLTData.lld'

#============= Step 1. Extract data in given domain ===================

echo 'Extracting data in defined domain'   

awk -v lonmin=$lonmin -v lonmax=$lonmax -v latmin=$latmin -v latmax=$latmax \
    '{ if ( $7 >= lonmin && $7 <= lonmax && $8 >= latmin && $8 <= latmax ) print }'\
    $Hermandata > $ExtractedDomain

#============ Step 2. Eliminate bad data and interpolated data =========

# In this step bad data are Eliminated 
# and a warning with the lines 
# that have bad data are printed in $BadaData file
# Interpolated data are also thrown away

echo 'Eliminating Bad and interpolated data'

awk -v GoodData=$GoodData -v BadData=$BadData -v shore=$shore\
  '{  if (NF != 15) {
        print "number of columns in line",NR,"is not equal to 14" > BadData
        print > BadData 
       }   
       else

       if ($13 !~ /MWL|MSL|LLWLT/) {
         print "line",NR,"does not have appropriate reference level"  > BadData
         print > BadData
       }         
       else 
       
       if ($10 < -10000) {
         print "check depth unit of line",NR > BadData
         print > BadData
       }       
       else

       if ($11!="GMT" && $11!="NETCDF" && $11!="NGDC" && $11!="IBCAO" && $14!="ABC" && $11!="CHSQUE" && $11!= "NSTOPO" && $14!="INTERP" ) {
         print > GoodData
       }
               
  }' $ExtractedDomain 
       
 
#======= Step 3. Separate data that should be corrected with WebTide ========

echo 'Separating Data to be corrected with WebTide'

awk -v MSLData=$MSLData -v LLWLTData=$LLWLTData \
   '{  if ($13 ~ /MWL|MSL/ || $10 < -200) {
         print $7,$8,$10 > MSLData
        }
        else
        if ($13 ~ /LLWLT/) {
          print $7,$8,$10 > LLWLTData
        }
     }' $GoodData
     



