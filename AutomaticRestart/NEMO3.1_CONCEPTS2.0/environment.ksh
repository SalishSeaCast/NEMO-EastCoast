#!/bin/ksh
set -ax
##	    environment.ksh

## ======================================================================
## ===============
## 2014/03/05 : JP Paquin
## Defines the environment and PATH for running NEMO with automatic 
## resubmission capacities for NEMO 3.1
## HEAVILY BASED ON MERCATOR VERSION
## ===============
##
## ======================================================================
## ===============
## Initialisations
## ===============
## MAXSUB        : automatic qsub up to MAXSUB lines in $CONFIG.db
## CONF_NAME     : Script used to run the simulation
## EXEC          : executable name
## NBPROC        : number of processors
## NDATEDEB      : date beginning of simulation
## *******  Begin Pegasus section  ******
## TIMELIMIT     : Time limit for Pegasus : DO NOT CHANGE
## MEMLIMIT      : Memory limit (Fred's)  : DO NOT CHANGE
## PEGASUS_QUEUE : Queue to use on Pegasus: DO NOT CHANGE
## SEAICE_MODULE : Sea ice module to use (LIM/CICE/NONE)
## PASDETEMPS    : Time step: 86400/PASDETEMPS must be integer
## OUTPUTFREQ    : Output frequency (in days!!!)
## ATMOS_DATA    : Atmospheric forcing (CGRF or ERAI)
## RST_FILE_FORMAT:Format of restart files (netcdf,dimg)
## OBC_ACTIVE    : Active OBC in the config
## ======================================================================
MAXSUB=37
CONF_NAME=CREG025_LIM
CONF_NBR=T17
CONFIG=${CONF_NAME}-${CONF_NBR}
SYSTEM_ABOVE=
NBPROC=144
#NBPROC_XIO=8     # XIO Server capacity not ready on Pegasus
NDATEDEB=20020102
# *******  Begin Pegasus section  ******
TIMELIMIT=10800   
MEMLIMIT=2000M
PEGASUS_QUEUE=dev
SEAICE_MODULE=LIM
PASDETEMPS=1080
OUTPUTFREQ=2
ATMOS_DATA=CGRF
RST_FILE_FORMAT=netcdf
OBC_ACTIVE='yes'

## ======================================================================
## SCRIPT OPTIONS  : initialisation and archivement
## 
## archive = Archive list.files built in POST-TREATEMENT.ksh 
## beg_rst = Begin from restart files (if first run or not N.B. Harcoded names RAPPAT.ksh)
## ======================================================================
archive=false
build_nc_mpp=false
beg_rst=true

## ======================================================================
## "MACHINECALCU" directories
## ======================================================================
## DIR_CALCU_EXE    : directory where is the executable			(rappat)
## DIR_CALCU_SRC    : directory where are model sources 		(rappat; ?USED?)
## DIR_CALCU_CTL    : directory where are control files of the run 	(namelists, *.xml not used here)
## DIR_CALCU_UTL    : directory where to find utilities
## ======================================================================
DIR_CALCU_EXE=/home/rpne/rpnejpp/MODEL_code/CONCEPTS_2.0.0/nemo3_1_cmc_${CONFIG}/modipsl/bin/${EC_ARCH}
DIR_CALCU_SRC=/home/rpne/rpnejpp/MODEL_code/CONCEPTS_2.0.0/nemo3_1_cmc_${CONFIG}/modipsl/modele/NEMO/WORK 
DIR_CALCU_CTL=/home/rpne/rpnejpp/MODEL_code/CONCEPTS_2.0.0/nemo3_1_cmc_${CONFIG}/CONTROL_FILES
DIR_CALCU_UTL=${HOME}/bin # JPP not used at the moment (reserve for post-processing)

## ======================================================================
## "MACHINEARCHI" directories
## ======================================================================
## DIR_ARCHI_FORC      : directory to stock the forcing fields		(rappat)
## DIR_ARCHI_INPUT     : directory for the input files			(rappat)
## DIR_ARCHI_CLIM      : directory for the climatology files 		(rappat)
## DIR_ARCHI_BIO       : directory for the input files for the biology
## DIR_ARCHI_RST_ABOVE : directory for restart of the system above
## DIR_ARCHI_OBC       : directory for the obc files			(rappat)
## DIR_ARCHI_OUTPUT    : directory to stock the results			(post_?)
## DIR_ARCHI_RST       : directory to stock the restarts		(post_?)
## ======================================================================
DIR_ARCHI_FORC=/fs/peg/data/rpne/rpnejpp/DATA/CGRF_at_Mercator/DATA
DIR_ARCHI_INPUT=/fs/peg/data/rpne/rpnejpp/INITIALISATION/CREG025
DIR_ARCHI_CLIM=/fs/peg/data/rpne/rpnejpp/INITIALISATION/CREG025
DIR_ARCHI_BIO=not_used
DIR_ARCHI_RST_ABOVE=not_used
DIR_ARCHI_OBC=/fs/peg/data/rpne/rpnejpp/INITIALISATION/CREG025/CREG025_IN_50lvls_FRED
CDIR=/fs/peg/data/rpne/rpnejpp
DIR_ARCHI_OUTPUT=$CDIR/${CONF_NAME}/${CONFIG}-S
DIR_ARCHI_RST=$CDIR/${CONF_NAME}/${CONFIG}-RST
DIR_ARCHI_CMB=$CDIR/${CONF_NAME}/${CONFIG}-OUT/CDF_COMB
DIR_ARCHI_ICE=/fs/peg/data/rpne/rpnejpp/INITIALISATION/CREG025

## ======================================================================
## "TEMP" directories
## ======================================================================
## DIR_TMP_HOME   : home temp directory 
## DIR_TMP_RUN    : directory to run
## DIR_TMP_RST    : directory for the restart
## DIR_TMP_OUT    : directory for the output
## DIR_TMP_DIMG   : directory for the dimg
## DIR_TMP_CDF    : directory for the cdf
## DIR_TMP_LOG    : directory for the log
## DIR_TMP_MOOR   : directory for the moorings
## DIR_TMP_DCT    : directory for the diadct
## DIR_TMP_INPUT  : directory for the input files
## DIR_TMP_FORC   : directory for the atmos forcing
## DIR_TMP_OBC    : directory for the obc files
## ======================================================================
DIR_TMP_HOME=/home/rpne/$USER/${CONF_NAME}
DIR_TMP_RUN=$DIR_TMP_HOME/DIRRUN_${CONFIG}
DIR_TMP_RST=$DIR_TMP_HOME/${CONFIG}-RST
DIR_TMP_OUT=$DIR_TMP_HOME/${CONFIG}-OUT
DIR_TMP_DIMG=$DIR_TMP_OUT/DIMG
DIR_TMP_MEAN=$DIR_TMP_HOME/${CONFIG}-MEAN
DIR_TMP_WORK=$DIR_TMP_OUT/WORK
DIR_TMP_CDF=$DIR_TMP_OUT/CDF
DIR_TMP_LOG=$DIR_TMP_OUT/LOG
DIR_TMP_MOOR=$DIR_TMP_OUT/MOORINGS
DIR_TMP_DCT=$DIR_TMP_OUT/DIADCT
DIR_TMP_INPUT=$DIR_TMP_HOME/${CONF_NAME}-IN
DIR_TMP_FORC=$DIR_TMP_RUN/ATMDATA
DIR_TMP_OBC=$DIR_TMP_HOME/OBCDATA

if [ ! -d $DIR_TMP_HOME ] ; then mkdir $DIR_TMP_HOME ; fi
if [ ! -d $DIR_TMP_RUN ] ; then mkdir $DIR_TMP_RUN ; fi
if [ ! -d $DIR_TMP_RST ] ; then mkdir $DIR_TMP_RST ; fi
if [ ! -d $DIR_TMP_OUT ] ; then mkdir $DIR_TMP_OUT ; fi
if [ ! -d $DIR_TMP_DIMG ] ; then mkdir $DIR_TMP_DIMG ; fi
if [ ! -d $DIR_TMP_MEAN ] ; then mkdir $DIR_TMP_MEAN ; fi
if [ ! -d $DIR_TMP_WORK ] ; then mkdir $DIR_TMP_WORK ; fi
if [ ! -d $DIR_TMP_CDF ] ; then mkdir $DIR_TMP_CDF ; fi
if [ ! -d $DIR_TMP_LOG ] ; then mkdir $DIR_TMP_LOG ; fi
if [ ! -d $DIR_TMP_MOOR ] ; then mkdir $DIR_TMP_MOOR ; fi
if [ ! -d $DIR_TMP_DCT ] ; then mkdir $DIR_TMP_DCT ; fi
if [ ! -d $DIR_TMP_INPUT ] ; then mkdir $DIR_TMP_INPUT ; fi
if [ ! -d $DIR_TMP_FORC ] ; then mkdir $DIR_TMP_FORC ; fi
if [ ! -d $DIR_TMP_OBC ] ; then mkdir $DIR_TMP_OBC ; fi

## ======================================================================
## FILES NAMES
## ======================================================================
EXEC_NAME=opa
BATMET_FIL=bathy_arctic025_mask.nc
COORD_FIL=coordinates.nc
MOOR_FIL=
SECT_FIL=
DCT_FIL=
TMX_K1=
TMX_M2=
TMX_MASK=
TMX_MASK=
SHLAT_FIL=
TEMP_FIL=CREG025_20020115.nc
SAL_FIL=CREG025_20020115.nc
U_FIL=CREG025_20020115.nc
V_FIL=CREG025_20020115.nc
SSH_FIL=CREG025_20020115.nc
RUNOFF_FIL=runoff_CREG025.nc 
CHLA_FIL=
INITICE_FIL=Ice_initialization_LIM_CONCEPTS.nc
BDY_FILU=
BDY_FILV=
BDY_FILT=
TIDE_FILT=
TIDE_FILU=
TIDE_FILV=
DUST_FIL=
RIVER_FIL=
obc_list='obcnorth.nc obcsouth.nc'  # Names Hardcoded in obcdta.F90 & obcdta_mer.F90
RESH_FIL=RPN_weights_025.nc
RESH_FIL_V=


## ======================================================================
## FUNCTIONS rappatrie, rappatrest and expatrie
## ======================================================================
	rappatrie ()
	{
	    if [ -f $2/$1 ] ; then
		ln -sf $2/$1 $1 ;
	    else
		cp $3 $1;
		cp $1 $2/$1 ;
	    fi ;
	}
	rappatrie2 ()
	{
	    if [ -f $3/$2 ] ; then
		echo  $3/$2 PRESENT;
		ln -sf $3/$2 $1 ;
	    else
		echo  $3/$2 NOT PRESENT ;
		tarfile=ERAinterim_$6_$5.tar ;
		cp $4/${tarfile}.gz $3 ;
		cd $3 ;
		gunzip $3/${tarfile}.gz ;
		tar xvf $3/${tarfile} ;
		cd ${DIR_TMP_RUN} ; 
		ln -sf  $3/$2 $1 ;
	    fi ;
	}
	rappatrest ()
	{
	    if [ -f $2 ] ; then
		ln -sf $2 $1 ;
	    else
		cp $3 $1 ;
	    fi ;
	}
        expatrie () {
            cp $2/$1 $3/$1 ;
	}
	expatrie_mv () {		 # JPP - move faster than copy
	    if [ ! -d $3 ] ; then
		echo 'mkdir directory: '$3
		mkdir -m 755 -p $3
	    fi
	    mv $2/$1 $3/$1 ;
	}

## ======================================================================
## FUNCTION j2d
## ======================================================================
	j2d ()
	{
#            set -axe
#	    typeset -i day month year tmpday centuries
#	    ((tmpday = $1 + 712164))
#	    ((centuries = (4 * tmpday - 1) / 146097))
#	    ((tmpday += (centuries - centuries/4)))
#	    ((year = (4 * tmpday - 1) / 1461))
#	    ((tmpday -= (1461 * year) / 4))
#	    ((month = (10 * tmpday - 5) / 306))
#	    ((day = tmpday - (306 * month + 5) / 10))
#	    ((month += 2))
#	    ((year += month/12))
#	    ((month = month % 12 + 1))
            set -axe
#           typeset -i day month year tmpday centuries
           let tmpday=$(($1+712164))
           let centuries=$(((4*$tmpday-1)/146097))
           let tmpday=$(($tmpday + ($centuries - $centuries/4)))
           let year=$(( (4 * $tmpday - 1) / 1461))
	   let tmpday=$(( $tmpday- (1461 * $year) / 4))
           let month=$(( (10 * $tmpday - 5) / 306))
           let day=$(( $tmpday - (306 * $month + 5) / 10))
           let month=$(( $month + 2))
           let year=$(( $year + $month/12))
           let month=$(( $month % 12 + 1))
	    echo ${year}`printf %2.2i ${month}``printf %2.2i ${day}`
	}

## ======================================================================
## FUNCTION d2j
## ======================================================================
	d2j ()
	{
	    typeset -i day month year tmpmonth tmpyear
	    ndate=$1
	    year=`echo ${ndate} | cut -c 1-4`
	    month=`echo ${ndate} | cut -c 5-6`
	    day=`echo ${ndate} | cut -c 7-8`
	    ((tmpmonth = 12 * year + month - 3))
	    ((tmpyear = tmpmonth / 12))
	    print $(((734*tmpmonth + 15)/24 - 2*tmpyear + tmpyear/4 - tmpyear/100 + tmpyear/400 + day - 712164))
	}

## ======================================================================
## Print informations relative to the run
## ======================================================================
echo $DIR_CALCU_CTL > ~/RUN_${CONF_NAME}-${CONF_NBR}.txt

echo
echo Configuration name : $CONFIG
echo Number of procs : $NBPROC
echo
if [ $beg_rst = "true" ] ; then
    echo Hot restart from the system $SYSTEM_ABOVE at the date $NDATEDEB
else
    echo Initialisation from climatology or analytical profiles
fi
echo
if [ $build_nc_mpp = "true" ] ; then
    echo The build_nc is multi proc
else
    echo The build_nc is mono proc
fi
echo
if [ $archive = "true" ] ; then
    echo The results will be archived
else
    echo The results will not be archived
fi
echo
echo home : $HOME
echo 
echo "MACHINECALCU" directories :
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo $DIR_CALCU_EXE    : directory where is the executable
echo $DIR_CALCU_SRC    : directory where are model sources
echo $DIR_CALCU_CTL    : directory where are control files of the run
echo $DIR_CALCU_UTL    : directory where to find utilities
echo
echo "MACHINEARCHI" directories :
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo $DIR_ARCHI_FORC      : directory to stock the forcing fields
echo $DIR_ARCHI_INPUT     : directory for the input files
echo $DIR_ARCHI_CLIM      : directory for the climatology files
echo $DIR_ARCHI_BIO       : directory for the input files for the biology
echo $DIR_ARCHI_RST_ABOVE : directory for restart of the system above
echo $DIR_ARCHI_OBC       : directory for the obc files
echo $DIR_ARCHI_OUTPUT    : directory to stock the results
echo $DIR_ARCHI_RST       : directory to stock the restarts
echo $DIR_ARCHI_CMB       : directory for combined netCDF results
echo $DIR_ARCHI_ICE       : directory for ice something ...
echo
echo "TEMP" directories :
echo ~~~~~~~~~~~~~~~~~~~~
echo $DIR_TMP_HOME   : home temp directory 
echo $DIR_TMP_RUN    : directory to run
echo $DIR_TMP_RST    : directory for the restart
echo $DIR_TMP_OUT    : directory for the output
echo $DIR_TMP_DIMG   : directory for the dimg
echo $DIR_TMP_CDF    : directory for the cdf
echo $DIR_TMP_LOG    : directory for the log
echo $DIR_TMP_MOOR   : directory for the moorings
echo $DIR_TMP_DCT    : directory for the diadct
echo $DIR_TMP_INPUT  : directory for the input files
echo $DIR_TMP_FORC   : directory for the atmos forcing
echo $DIR_TMP_OBC    : directory for the obc files
echo
echo Run features :
echo ~~~~~~~~~~~~~~
echo Executable name : $EXEC_NAME
echo Bathymetry : $BATMET_FIL
echo Coordinates : $COORD_FIL
echo Position.moor : $MOOR_FIL
echo Position.sect : $SECT_FIL
echo Section for diadct : $DCT_FIL
echo K1 tidal mixing : $TMX_K1
echo M2 tidal mixing : $TMX_M2
echo Mask tidal mixing : $TMX_MASK
echo Shlat : $SHLAT_FIL
echo Clim T : $TEMP_FIL
echo Clim U : $SAL_FIL
echo Runoff : $RUNOFF_FIL
echo Chlorophyl : $CHLA_FIL
echo Bdy U : $BDY_FILU
echo Bdy V : $BDY_FILV
echo Bdy T : $BDY_FILT
echo Tide T : $TIDE_FILT
echo Tide U : $TIDE_FILU
echo Tide V : $TIDE_FILV
echo List of the open boundaries : $obc_list
echo
