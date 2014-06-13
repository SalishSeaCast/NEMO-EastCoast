#!/bin/ksh
## ======================================================================
##            SCRIPT TO SUCCESSIVELY SUBMIT JOBS ON BALTIC
##	      RAPPAT.ksh
## ======================================================================
echo "\n############## ${Job_backend} CONFIGURATION ###########\n"
. ~/.profile_usr
echo "\n############## ${Job_backend} CONFIGURATION (end) ###########\n"

# JPP : Text file ~/RUN_${CONF_NAME}-${CONF_NBR}.txt was created by encivronment.ksh
#       to store de directory where the control scripts are located
#       because following ord_soumet, the directory is back to be $HOME
#       It is patchy, requires manual modification but it seems to work !!!
file_RUN=~/RUN_CREG025_LIM-T17.txt 
if [ -f ${file_RUN} ] ; then 
  direnv=`tail -1 ${file_RUN}`
  cd $direnv
fi

set -ex 
. ./environment.ksh

## ======================================================================
## PREPARE THE RUN 
## ======================================================================
cd $DIR_TMP_RUN
rm -rf $DIR_TMP_RUN/*	# Clean DIR_TMP_RUN in automated sequence

## ======================================================================
## Copy the executable 
## ======================================================================
if [ ! -f ${DIR_CALCU_EXE}/$EXEC_NAME ] ; then
    echo ${DIR_CALCU_EXE}/$EXEC_NAME 'IS MISSING !!!!!'
    exit
else
#  JPP - XIO Server not available on Pegasus
#       if [ ${NBPROC_XIO} -ge 0 ] ; then
#          cp ${DIR_CALCU_CTL}/iodef.xml .
#          cp ${DIR_CALCU_CTL}/field_def.xml .
#          cp ${DIR_CALCU_CTL}/domain_def.xml .
#          cp ${XIOS_BIN}/xios_server.exe .
#          #cp ${DIR_CALCU_CTL}/xios_server.exe .
#       fi

    cp ${DIR_CALCU_EXE}/$EXEC_NAME ./opa
    chmod 777 ./opa
fi

## ======================================================================
## Copy the namelist
## ======================================================================
if [ ! -f $DIR_CALCU_CTL/namelist-${CONFIG} ] ; then
    echo $DIR_CALCU_CTL/namelist-${CONFIG} 'IS MISSING !!!!!' 
    exit
else
    cp $DIR_CALCU_CTL/namelist-${CONFIG} namelist
fi

# Get Sea Ice Namelist
# LIM2
if [ $SEAICE_MODULE == 'LIM' ] ; then
    if [ ! -f $DIR_CALCU_CTL/namelist_lim-${CONFIG} ] ; then
      echo $DIR_CALCU_CTL/namelist_lim-${CONFIG} 'IS MISSING !!!!!'; exit
    else
      cp $DIR_CALCU_CTL/namelist_lim-${CONFIG} namelist_ice
    fi
fi
# CICE
if [ $SEAICE_MODULE == 'CICE' ] ; then
    if [ ! -f $DIR_CALCU_CTL/namelist_cice-${CONFIG} ] ; then
      echo $DIR_CALCU_CTL/namelist_cice-${CONFIG} 'IS MISSING !!!!!'; exit
    else
      cp $DIR_CALCU_CTL/namelist_cice-${CONFIG} namelist_ice
    fi
fi


## ======================================================================
## Copy the $CONFIG.db
## ======================================================================
if [ ! -f $DIR_CALCU_CTL/$CONFIG.db ] ; then
    echo $DIR_CALCU_CTL/$CONFIG.db  'IS MISSING !!!!!' 
    exit
else
    cp $DIR_CALCU_CTL/$CONFIG.db .
fi


## ======================================================================
## Date and flag
## Modify the $CONFIG.db for automatic re-submitting 
## ======================================================================
line=` tail -1 $CONFIG.db `
no=`echo $line|cut -d' ' -f1`
nit000=` echo $line|cut -d' ' -f2`
nitend=` echo $line|cut -d' ' -f3`
if [ $no -ne 1 ] ; then
    line2=` tail -2 $CONFIG.db`
    ndaterest=` echo $line2|cut -d' ' -f4`
    ndate0=$(/bin/date "+%Y%m%d" --date="${ndaterest} 1 days")
    # JPP restnum : time step of the restart...
    clstep=` echo $line2|cut -d' ' -f3`
    typeset -Z8 clstep=$clstep
    flag_rst=.true.
    flag_ramp=.false.
    flag_rst_trc=.true.
    rsttype=2    
else
    ndate0=$NDATEDEB
    if [ $beg_rst == "true" ]; then
        flag_rst=.false.
        flag_ramp=.true.
        flag_rst_trc=.false.
    else
        flag_rst=.false.
        flag_ramp=.true.
        flag_rst_trc=.false.
    fi
    rsttype=0
fi


## ======================================================================
## Time step and output frequencies
## ======================================================================
freqo=` echo 1 | awk "{ a=int( $OUTPUTFREQ *86400 / $PASDETEMPS) ; print a }"`

## ======================================================================
## Calendar
## ======================================================================
ndays=` echo 1 | awk "{ a=int( ($nitend - $nit000 +1)*$PASDETEMPS /86400.) ; print a }"`
ydeb=`date --date="${ndate0}" +%Y`
mdeb=`date --date="${ndate0}" +%m`
ddeb=`date --date="${ndate0}" +%d `
ndastpfin=$(/bin/date "+%Y%m%d" --date="${ndate0} ${ndays} days")
yfin=`date --date="${ndastpfin}" +%Y`
mfin=`date --date="${ndastpfin}" +%m`
dfin=`date --date="${ndastpfin}" +%d `
julstart=$(d2j ${ndate0})
julstop=$(d2j ${ndastpfin})
echo $ndays days to run, starting $ndate0 "("$julstart")" ending $ndastpfin "("$julstop")"


## ======================================================================
## Rappat : Copy the files adapted to the config
## ======================================================================
## ======================================================================
## Bathymetry
## ======================================================================
nomini=$BATMET_FIL
nom=bathy_meter.nc
rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi

## ======================================================================
## Coordinate
## ======================================================================
nomini=$COORD_FIL
nom=coordinates.nc
rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi

## ======================================================================
## Climatology
## ======================================================================
clim98=9999
if [ ${clim98} -ne 9999 ] ; then  ## OPTION NOT USED IF 9999
   if [ ${clim98} -eq 1 ] ; then
	# Temperature
        ((num=mdeb-1))
        if [  $num -eq 0 ] ; then num=12  ; fi
        typeset -Z2 num=${num}
        nomini=Levitus_p2.1_Tpot_CREG025_L75m${num}.nc
        nom=Levitus_p2.1_Tpot_CREG025_L75_m${num}.nc
        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_CLIM/$nomini
        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi

        ((num=mdeb))
        typeset -Z2 num=${num}
        nomini=Levitus_p2.1_Tpot_CREG025_L75m${num}.nc
        nom=Levitus_p2.1_Tpot_CREG025_L75_m${num}.nc
        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_CLIM/$nomini
        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi

        ((num=mdeb+1))
        if [  $num -eq 13 ] ; then num=1  ; fi
        typeset -Z2 num=${num}
        nomini=Levitus_p2.1_Tpot_CREG025_L75m${num}.nc
        nom=Levitus_p2.1_Tpot_CREG025_L75_m${num}.nc
        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_CLIM/$nomini
        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi

        # Salinity
        ((num=mdeb-1))
        if [  $num -eq 0 ] ; then num=12  ; fi
        typeset -Z2 num=${num}
        nomini=Levitus_p2.1_S_correc_CREG025_L75m${num}.nc
        nom=Levitus_p2.1_S_correc_CREG025_L75_m${num}.nc
        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_CLIM/$nomini
        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi

        ((num=mdeb))
        typeset -Z2 num=${num}
        nomini=Levitus_p2.1_S_correc_CREG025_L75m${num}.nc
        nom=Levitus_p2.1_S_correc_CREG025_L75_m${num}.nc
        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_CLIM/$nomini
        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi

        ((num=mdeb+1))
        if [  $num -eq 13 ] ; then num=1  ; fi
        typeset -Z2 num=${num}
        nomini=Levitus_p2.1_S_correc_CREG025_L75m${num}.nc
        nom=Levitus_p2.1_S_correc_CREG025_L75_m${num}.nc
        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_CLIM/$nomini
        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
   else
	# Temperature
        nomini=Levitus05_potemp05_1m_z75_nomask_smoothed.nc
        nom=data_1m_potential_temperature_nomask.nc
        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_CLIM/$nomini
        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
       
 	# Salinity
        nomini=Levitus05_salin05_1m_z75_nomask_smoothed.nc
        nom=data_1m_salinity_nomask.nc
        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_CLIM/$nomini
        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
   fi
fi

## ======================================================================
## Ice initial conditions
## ======================================================================
#inice=`grep ln_limini namelist_ice | tail -1 | awk -F" " '{print $3}'`
#if [ ${inice} == ".true." ] ; then
if [ $no -eq 1 ] ; then
  if   [ $SEAICE_MODULE == 'LIM' ] ; then 
      nomini=$INITICE_FIL
      nom='Ice_initialization.nc'
      rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_ICE/$nomini
      if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
  elif [ $SEAICE_MODULE == 'CICE' ] ; then
      nomini=$INITICE_FIL
      nom='init_cice.nc'
      rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_ICE/$nomini
      if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
  else
      echo 'NO SEA ICE MODULE'
  fi
fi

## ======================================================================
## Chlorophyl 					JPP : NOT USED
## ======================================================================
#nchla=`grep nn_chldta  namelist | tail -1 | awk -F" " '{print $3}'`
#if [ $nchla -ge 1 ] ; then
#    nomini=$CHLA_FIL
#    #nom='kpar.nc'
#    nom='chlorophyll.nc'
#    rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
#    if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
#fi
## ======================================================================
## Runoff (precipitation form) 			JPP : NOT USED (see below)
## ======================================================================
#nrunof=`grep nrunoff  namelist | tail -1 | awk -F" " '{print $3}'`
#nomini=$RUNOFF_FIL
#nom=runoff_1m_nomask.nc
#if [ $nrunof -ge 1 ] ; then
#    rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
#    if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
#fi
## ======================================================================
## Runoff (bdy form)                    	JPP : NOT USED (see below)
## ======================================================================
#grep 'key_bdy' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#if [ $? -eq 0 ] ; then
#    nomini=$BDY_FILT
#    nom='bdydta_T.nc'
#    rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
#    nomini=$BDY_FILU
#    nom='bdydta_U.nc'
#    rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
#    nomini=$BDY_FILV
#    nom='bdydta_V.nc'
#    rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
#fi
## ======================================================================
## Runoff
## ======================================================================
nomini=$RUNOFF_FIL
nom='runoff_1m_nomask.nc'
rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi

## ======================================================================
## physics obc files 
## ======================================================================
## PEGASUS
# Hardcoded names for obc files obc${direction}_[TS,U,V].nc
# in obcdta.F90 & obcdta_mer.F90
if [ $OBC_ACTIVE = 'yes' ] ; then
#    ln -s $DIR_ARCHI_OBC .
    for file in ${obc_list} ; do
      fileobc=`echo $file | sed 's/.nc//g'`
      rappatrie ${fileobc}.nc $DIR_TMP_OBC $DIR_ARCHI_OBC/${fileobc}.nc
      if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
      ln -s ${fileobc}.nc ${fileobc}_TS.nc
      ln -s ${fileobc}.nc ${fileobc}_U.nc
      ln -s ${fileobc}.nc ${fileobc}_V.nc
    done
    ## ==================================================================
    ## physics obc files - tidal part JPP NOT USED
    ## ==================================================================
    #grep 'key_tide' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
    #if [ $? -eq 0 ] ; then
    #    ln -s ${DIR_ARCHI_INPUT}/tidal_forcing .
    #fi
    #ln -s ${DIR_ARCHI_INPUT}/rnf_forcing .
    #ln -s ${DIR_ARCHI_INPUT}/coordinates.bdy.nc .
fi
## MERCATOR
#grep 'key_bdy' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#if [ $? -eq 0 ] ; then
#     ln -s $DIR_ARCHI_OBC .
#    ## ==================================================================
#    ## physics obc files - tidal part
#    ## ==================================================================
#    grep 'key_tide' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#    if [ $? -eq 0 ] ; then
#        ln -s ${DIR_ARCHI_INPUT}/tidal_forcing .
#    fi
#    ln -s ${DIR_ARCHI_INPUT}/rnf_forcing .
#    ln -s ${DIR_ARCHI_INPUT}/coordinates.bdy.nc .
#fi
## ======================================================================
## Restart
## ======================================================================
if [ $no -eq 1 ] ; then
    if [ $beg_rst == "false" ] ; then
      echo "Initialization from climatology or an analytical profile"
      echo "JPP: Do not think this works !!!"  ; exit
    else
# PEGASUS
        # Temperature
        nomini=$TEMP_FIL
        nom=IC_T.nc
        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_CLIM/$nomini
        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
        # Salinity
        nomini=$SAL_FIL
        nom=IC_S.nc
        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_CLIM/$nomini
        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
        # SSH
        nomini=$SSH_FIL
        nom=IC_SSH.nc
        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_CLIM/$nomini
        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
        # U velocity
        nomini=$U_FIL
        nom=IC_U.nc
        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_CLIM/$nomini
        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
        # V Velocity
        nomini=$V_FIL
        nom=IC_V.nc
        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_CLIM/$nomini
        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
    fi
else
   # JPP : Restart file names (or prefix) are defined in 
   #       namelist (cn_ocerst_in) and namelist_ice (cn_icerst_in) 
   #       and must fit the prefix here! 
   jproc=0
   NBPROCM1=`expr $NBPROC - 1`
   while [ $jproc -le $NBPROCM1 ] ; do
      if [ $RST_FILE_FORMAT == 'netcdf' ] ; then 
      # ------   RESTART IN NETCDF FORMAT   -----
      #  nb of Proc from 0 to NBPROC-1 
      #  jproc must be 4 digits
        typeset -Z4 jproc

	# OCEAN
        nomrest=restart_${jproc}.nc
        nomrestart=${CONFIG}_${clstep}_restart_${jproc}.nc
  	ln -s $DIR_ARCHI_RST/$nomrestart $nomrest
	# SEA ICE 
        nomrest=restart_ice_${jproc}.nc
        nomrestart=${CONFIG}_${clstep}_restart_ice_${jproc}.nc
	ln -s $DIR_ARCHI_RST/$nomrestart $nomrest

      elif [ $RST_FILE_FORMAT == 'dimg' ] ; then 
      # -----   RESTART IN DIMG FORMAT (requires key_dimgout)   -----
      #   nb of Proc from 1 to NBPROC
      #   jproc must be 5 digits
      #   !!!UNTESTED AT THE MOMENT!!!
        typeset -Z5 jproc
        typeset -Z5 jprocm1
        nomrest=restart_${jprocm1}.dimg
        nomrestart=${CONFIG}_${clstep}_restart_${jprocm1}.nc
        rappatrest $nomrest $DIR_TMP_RST/$nomrestart $DIR_ARCHI_RST/$nomrestart
        if [ $? -ne 0 ] ; then echo "$nomrestart non trouve" ; exit ; fi

        nomrest=restart_ice_${jprocm1}.nc
        nomrestart=${CONFIG}_${clstep}_restart_ice_${jprocm1}.nc
        rappatrest $nomrest $DIR_TMP_RST/$nomrestart $DIR_ARCHI_RST/$nomrestart
        if [ $? -ne 0 ] ; then echo "$nomrestart non trouve" ; exit ; fi
      fi
      let jproc=$jproc+1
    done
fi
 
## ======================================================================
## Restart OBC
## ======================================================================
echo 'JPP debug Restart OBC no='${no}
if [ $no -ne 1 ] ; then
  if [ ${OBC_ACTIVE} = 'no' ] ; then
    nom=restart.obc.output
    nomini=${CONFIG}_${clstep}_restart_obc.output
    ln -s $DIR_ARCHI_RST/$nomini $nom
    if [ $? -ne 0 ] ; then echo "$nomrestart non trouve" ; exit ; fi
  elif [ ${OBC_ACTIVE} = 'yes' ] ; then
    nom=restart.obc.output
    nomini=${CONFIG}_${clstep}_restart_obc.output
    ln -s $DIR_ARCHI_RST/$nomini $nom
  fi
fi


## ======================================================================
## Atmospherical forcings
## ======================================================================
if   [ $ATMOS_DATA == 'ERAI' ] ; then  
> atm.list # Check atm. variables needed for the run
  echo BULKHUMI >> atm.list
  echo BULKTAIR >> atm.list
  echo FLCOPRE_v2  >> atm.list
  echo BULKU10M >> atm.list
  echo BULKV10M >> atm.list
  echo FLCOSSRD_v2 >> atm.list
  echo FLCOSTRD_v2 >> atm.list
  #echo PRES >> atm.list
  echo SF >> atm.list
elif [ $ATMOS_DATA == 'CGRF' ] ; then 
> atm.list # Check atm. variables needed for the run
  echo q2 >> atm.list
  echo t2 >> atm.list
  echo precip  >> atm.list
  echo u10 >> atm.list
  echo v10 >> atm.list
  echo qsw >> atm.list
  echo qlw >> atm.list
  #echo PRES >> atm.list
  echo snow >> atm.list
else
  echo $ATMOS_DATA ' NOT KNOWN... EXIT' ; exit
fi
list_forc=`awk -F" " '{printf "%s " ,$1}' atm.list`
echo $list_forc
#rm -f atm.list

if [ ! -d $DIR_TMP_FORC ] ; then
  mkdir -p $DIR_TMP_FORC
fi


for typeforc in $list_forc ; do  # Loop on input files
    case $typeforc in
      # ERAI
      BULKHUMI | BULKTAIR | BULKU10M | BULKV10M  ) freq=3h ;;
      FLUXPRE | FLCOPRE_v2 | SF | FLUXSSRD | FLUXSTRD | FLCOSSRD_v2 | FLCOSTRD_v2 | PRES   ) freq=24h ;;
      # CGRF 
      q2 | t2 | u10 | v10       ) freq=3h  ;;
      precip | snow | qsw | qlw ) freq=24h ;;
    * ) echo 'UNKNOWN ATMOSPHERIC VARIABLE... EXIT' ; exit ;;
    esac

    # Files/variables names
    typeforc1=$typeforc
    if [ $ATMOS_DATA == 'ERAI' ] ; then 
      if [ $typeforc = 'FLCOPRE_v2' ] ; then typeforc1='FLCOPRE' ; fi 
      if [ $typeforc = 'FLCOSSRD_v2' ] ; then typeforc1='FLCOSSRD' ; fi 
      if [ $typeforc = 'FLCOSTRD_v2' ] ; then typeforc1='FLCOSTRD' ; fi 
      prefix=ERAinterim_${typeforc1}
    elif [ $ATMOS_DATA == 'CGRF' ] ; then
      prefix=$typeforc
    fi

    ((jul=julstart-1))
 
    # --- Basta, on y va par la fonction date... fini l'ostie de niaisage!!
    if [ $no -eq 1 ] ; then
      datec=${ndate0}
    else                     
      #-- Add link to previous day if restart (might avoid warnin msg)
      datec=$(/bin/date "+%Y%m%d" --date="${ndate0} -1 days")
    fi

    while [ $datec -le $ndastpfin ] ; do 
      year=`date  --date="${datec}" +%Y`
      month=`date --date="${datec}" +%m`
      day=`date   --date="${datec}" +%d `

      datec=$(/bin/date "+%Y%m%d" --date="${datec} 1 days")
    
      nomfic=${prefix}_${year}${month}${day}.nc
      nomfic2=${prefix}_y${year}m${month}d${day}.nc

      if [ $ATMOS_DATA == 'CGRF' ] ; then   # CGRF file format is ${var}_yYYYYmMMdDD.nc
        nomfic=${nomfic2}     
      fi 

      if [ ! -f $DIR_ARCHI_FORC/$typeforc/${freq}/${nomfic} ] ; then 
         echo $DIR_ARCHI_FORC/$typeforc/${freq}/${nomfic} ' DOES NOT EXIST ... EXIT' ; exit
      else
        ln -s $DIR_ARCHI_FORC/$typeforc/${freq}/${nomfic} ${DIR_TMP_FORC}/${nomfic2}
      fi
    done
done
## ======================================================================
## reshape - Online interpolation of Atm forcing
## ======================================================================
if [ $ATMOS_DATA == 'CGRF' ] ; then
  ## Only 1 file required for CGRF (scalar and vectorial)
  nomini=$RESH_FIL
  nom=$RESH_FIL
  rappatrie $nom $DIR_TMP_FORC $DIR_ARCHI_INPUT/$nomini
  if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi

elif   [ $ATMOS_DATA == 'ERAI' ] ; then
  ## Scalar fields 
  nomini=$RESH_FIL
  nom=reshape_ERAi_CREG025_bilin.nc
  rappatrie $nom $DIR_TMP_FORC $DIR_ARCHI_INPUT/$nomini
  if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
  ## Vectorial fields
  nomini=$RESH_FIL_V
  nom=reshape_ERAi_CREG025_bicubic.nc
  rappatrie $nom $DIR_TMP_FORC $DIR_ARCHI_INPUT/$nomini
  if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
fi

## **********************************************************************
## **********************************************************************
## **********  BEGINNING OF SECTION UNUSED FOR CONCEPTS CODE   **********
## ======================================================================
## Moorings : JPP NOT USED
## ======================================================================
#moor=`grep ln_diamoor namelist | tail -1 | awk -F" " '{print $3}'`
#if [ ${moor} == ".true." ] ; then
#    nomini=$MOOR_FIL
#    nom='position.moor'
#    rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
#    if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
#fi
## ======================================================================
## ahmcoef JPP : NOT USED 
## ======================================================================
#nomini=ahmcoef
#nom=ahmcoef
#rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
#if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
## ======================================================================
## Diadct JPP : NOT USED
## ======================================================================
#grep 'key_diadct' $DIR_CALCU_SRC/cpp_${CONFIG}.fcm
#if [ $? -eq 0 ] ; then
#    nomini=$DCT_FIL
#    nom='section_ijglobal.diadct'
#    rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
#    if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
#fi
## ======================================================================
## shlat JPP : NOT USED
## ======================================================================
#shlat=`grep shlat namelist | tail -1 | awk -F" " '{print $3}'`
#if [ $shlat = "-100" ] ; then
#    nomini=$SHLAT_FIL
#    nom='shlat2d.nc'
#    rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
#    if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
#fi
## ======================================================================
## Tidal mixing files JPP - NOT USED
## ======================================================================
#grep 'key_zdftmx' $DIR_CALCU_SRC/cpp_${CONFIG}.fcm
#if [ $? -eq 0 ] ; then
#    nomini=$TMX_K1
#    nom='K1rowdrg.nc'
#    rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
#    if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
#    nomini=$TMX_M2
#    nom='M2rowdrg.nc'
#    rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
#    if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
#    nomini=$TMX_MASK
#    nom='mask_itf.nc'
#    rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_INPUT/$nomini
#    if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
#fi
## ======================================================================
## PISCES Climatology                           JPP : NOT USED
## ======================================================================
#grep 'key_top' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#kcount=$?
#grep 'key_pisces' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#((kcount=${kcount}+$?))
#grep 'key_dtatrc' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#((kcount=${kcount}+$?))
#if [ ${kcount} -eq 0 ] ; then
#    for i in Alkalini BFe BSi CaCO3 DCHL DFe DIC DOC DSi Fer GOC NCHL NFe NH4 NO3 O2 PHY2 PHY PO4 POC SFe Si ZOO2 ZOO ; do
#        nomini=ORCA1_${i}_${CONF_NAME}.nc
#        nom=data_1m_${i}_nomask.nc
#        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_BIO/$nomini
#        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
#    done
#fi
## ======================================================================
## PISCES obc files                             JPP : NOT USED
## ======================================================================
#grep 'key_top' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#kcount=$?
#grep 'key_pisces' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#((kcount=${kcount}+$?))
#grep 'key_obc_mer' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#((kcount=${kcount}+$?))
#if [ ${kcount} -eq 0 ] ; then
#    obctype=$obc_list
#    for OBC in ${obc_list} ; do
#        nomini=obc${OBC}_BIO.nc
#        nom=obc${OBC}_BIO.nc
#        rappatrie $nom $DIR_TMP_OBC $DIR_ARCHI_BIO/$nomini
#        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
#    done
#fi
## ======================================================================
## PISCES static files JPP - NOT USED
## ======================================================================
#grep 'key_top' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#kcount=$?
#grep 'key_pisces' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#((kcount=${kcount}+$?))
#if [ ${kcount} -eq 0 ] ; then
#    ln_dust=`grep ln_dustfer  namelist_pisces | tail -1 | awk -F" " '{print $3}'`
#    nomini=$DUST_FIL
#    nom=dust.orca.nc
#    if [ $ln_dust == ".true." ] ; then
#        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_BIO/$nomini
#        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
#    fi
#    ln_riv=`grep ln_river  namelist_pisces | tail -1 | awk -F" " '{print $3}'`
#    nomini=$RIVER_FIL
#    nom=river.orca.nc
#    if [ $ln_riv == ".true." ] ; then
#        rappatrie $nom $DIR_TMP_INPUT $DIR_ARCHI_BIO/$nomini
#        if [ $? -ne 0 ] ; then echo "$nomini non trouve" ; exit ; fi
#    fi
#fi
## **********      END OF SECTION UNUSED FOR CONCEPTS CODE     **********
## **********************************************************************
## **********************************************************************


## ======================================================================
## namelist update
## ======================================================================
sed -e "s/NUMERO_DE_RUN/$no/"        \
    -e "s/EXPER/\"$CONFIG\"/"        \
    -e "s/NIT000/$nit000/"           \
    -e "s/NITEND/$nitend/"           \
    -e "s/NDATE0/$ndate0/"           \
    -e "s/BREST/$flag_rst/"          \
    -e "s/RAMP/$flag_ramp/"          \
    -e "s/DATERST/$rsttype/"         \
    -e "s/FREQO/$freqo/"             \
    -e "s/PASDETEMPS/$PASDETEMPS/"   \
    -e "s%PATH_RST%${DIR_TMP_RST}%"  \
    -e "s%PATH_OUT%${DIR_TMP_DIMG}%" \
    -e "s%PATH_MOOR%${DIR_TMP_MOOR}%" namelist > namelist1
mv namelist1 namelist

## ======================================================================
## Soumission
## ======================================================================
# Text file containing the $DIR_CALCU_CTL for the $CONF_NAME.ksh to read in
echo $DIR_CALCU_CTL > ${file_RUN}

cd $DIR_CALCU_CTL

ord_soumet $CONF_NAME.ksh -mach pegasus -cpus ${NBPROC} -mpi -cm ${MEMLIMIT} -t ${TIMELIMIT} \
                          -listing ${DIR_CALCU_CTL} -q ${PEGASUS_QUEUE}

## ======================================================================
## End of the rappat script
## ======================================================================
