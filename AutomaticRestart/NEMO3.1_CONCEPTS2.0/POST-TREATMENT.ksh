#!/bin/ksh
## ======================================================================
##            SCRIPT TO SUCCESSIVELY SUBMIT JOBS ON BALTIC
##            POST-TREATMENT.ksh
## ======================================================================
# JPP : Text file ~/RUN_${CONF_NAME}-${CONF_NBR}.txt was created by encivronment.ksh
#       to store de directory where the control scripts are located
#       because following ord_soumet, the directory is back to be $HOME
#       It is patchy, requires manual modification but it seems to work !!!
file_RUN=~/RUN_CREG025_LIM-T17.txt
if [ -f ${file_RUN} ] ; then                        
  direnv=`tail -1 ${file_RUN}`
  cd $direnv
fi

echo "\n############## ${Job_backend} CONFIGURATION ###########\n"
. ~/.profile_usr
echo "\n############## ${Job_backend} CONFIGURATION (end) ###########\n"

set -ex
. ./environment.ksh

## ======================================================================
## $CONFIG.db update
## ======================================================================
cd $DIR_TMP_LOG

line=` tail -1 $CONFIG.db `
no=`echo $line|cut -d' ' -f1`
nit000=` echo $line|cut -d' ' -f2`
nitend=` echo $line|cut -d' ' -f3`
typeset -Z8 clstep=${nitend}
if [ $no -ne 1 ] ; then
    line2=` tail -2 $CONFIG.db`
    ndaterest=` echo $line2|cut -d' ' -f4`
    ndate0=$(/bin/date "+%Y%m%d" --date="${ndaterest} 1 days")
else
    ndate0=$NDATEDEB
fi
#rdt=` grep rn_rdt namelist | tr -d ' ' | cut -c1-20 | grep rn_rdt= | cut -c8-50 | cut -d, -f1 | cut -d. -f1 `
rdt=$PASDETEMPS
let nbiter=`expr 86400/$rdt`
let nbiter_hf=`expr ${nbiter}/24`
ndays=` echo 1 | awk "{ a=int( ($nitend - $nit000 +1)*$rdt /86400.) ; print a }"`
ndastpfin=$(/bin/date "+%Y%m%d" --date="${ndate0} ${ndays} days")
ndatedm1=$(/bin/date "+%Y%m%d" --date="${ndastpfin} -1 days")

let no=` expr $no+1`
aammdd=${ndatedm1}
last1=`tail -1 $CONFIG.db`
ncol=`echo $last1|wc -w`

if [ $ncol = 3 ] ; then
    sed -e "s/$last1/$last1\ $aammdd/" $CONFIG.db > tmpdb
    mv -f tmpdb $CONFIG.db
else
    echo fichier db deja a jour de la date $ncol
fi
let dif=` expr $nitend-$nit000+1`
let nit000=` expr $nitend+1`
let nitend=` expr $nitend+dif`
echo $no $nit000 $nitend >> $CONFIG.db
mv $CONFIG.db $DIR_CALCU_CTL/.

> $DIR_TMP_OUT/list.files

## ======================================================================
## Rename log files
## ======================================================================
#  Keep a copy of the namelist in $DIR_CALCU_CTL...
no_ori=$no
let no_ori=` expr $no_ori-1`
cp namelist     ${DIR_CALCU_CTL}/namelist_${CONFIG}_${no_ori}
if [ $SEAICE_MODULE == 'LIM' ] ; then
  cp namelist_ice $DIR_CALCU_CTL/namelist_lim-${CONFIG}_${no_ori} 
fi
# CICE
if [ $SEAICE_MODULE == 'CICE' ] ; then
  cp namelist_ice $DIR_CALCU_CTL/namelist_cice-${CONFIG}_${no_ori} 
fi


nomenclature_run=${ndate0}_${ndastpfin}_R${ndate0}
for i in namelist* ocean.output* time.step solver.stat ; do
    mv ${i} ${CONFIG}_${i}_${nomenclature_run}
done
nomtar=tarfile.${CONFIG}_${nomenclature_run}_annex
tar cvf $nomtar *${nomenclature_run}
echo ${nomtar} $DIR_TMP_LOG $DIR_ARCHI_OUTPUT >> $DIR_TMP_OUT/list.files

## ======================================================================
## Diadct
## ======================================================================
#grep '$(prefix)key_diadct' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#if [ $? -eq 0 ] ; then
#    name=${CONFIG}_*_diadct
#    nomdct=${CONFIG}_DIADCT_${nomenclature_run}
#    mv ${DIR_TMP_DIMG}/$name ${DIR_TMP_DCT}/$nomdct
#    echo ${nomdct} $DIR_TMP_DCT $DIR_ARCHI_OUTPUT >> $DIR_TMP_OUT/list.files
#fi

## ======================================================================
## Diamoor
## ======================================================================
#moor=`grep ln_diamoor ${DIR_CALCU_CTL}/namelist-${CONFIG} | tail -1 | awk -F" " '{print $3}'`
#if [ ${moor} == ".true." ] ; then
#    cd $DIR_TMP_MOOR
#    cp $DIR_CALCU_UTL/build_nc .
#    cp $DIR_CALCU_CTL/namelistio_standard namelistio
#    sed -e "s/CONFIG/$CONFIG/"    \
#        -e "s/RDT/$rdt/"          \
#        -e "s/NDATE0/$ndate0/"    \
#        -e "s/FREQ/$nbiter/"      \
#        -e "s/FQHF/$nbiter_hf/"   \
#        -e "s/REC3D/false/"       namelistio > namelist1
#    mv namelist1 namelistio
#    ln -s $DIR_TMP_INPUT/position.moor .
#    ./build_nc 'moorings'
#    nomtar=${CONFIG}_MOORINGS_${nomenclature_run}
#    for name in *.mooring.nc ; do
#        mv ${name} ${CONFIG}_${nomenclature_run}_${name}
#    done
#    tar cvf $nomtar ${CONFIG}_${nomenclature_run}*.nc
#    echo ${nomtar} $DIR_TMP_MOOR $DIR_ARCHI_OUTPUT >> $DIR_TMP_OUT/list.files
#    rm -f *.mooring* build_nc namelistio position.moor build_nc_moorings_OK
#fi
#
## ======================================================================
## Restarts
## ======================================================================
jproc=0
NBPROCM1=`expr $NBPROC - 1`
while [ $jproc -le $NBPROCM1 ] ; do
  if [ $RST_FILE_FORMAT == 'netcdf' ] ; then
  # ------   RESTART IN NETCDF FORMAT   -----
  #  nb of Proc from 0 to NBPROC-1 
  #  jproc must be 4 digits
    typeset -Z4 jproc
    nomrestart=${CONFIG}_${clstep}_restart_${jproc}.nc
    expatrie_mv $nomrestart $DIR_TMP_RST $DIR_ARCHI_RST 
    if [ $? -ne 0 ] ; then echo "$nomrestart non trouve" ; exit ; fi

    nomrestart=${CONFIG}_${clstep}_restart_ice_${jproc}.nc
    expatrie_mv $nomrestart $DIR_TMP_RST $DIR_ARCHI_RST
    if [ $? -ne 0 ] ; then echo "$nomrestart non trouve" ; exit ; fi

  elif [ $RST_FILE_FORMAT == 'dimg' ] ; then
  # -----   RESTART IN DIMG FORMAT (requires keydimgout)   -----
  #   nb of Proc from 1 to NBPROC
  #   jproc must be 5 digits
  echo !!!UNTESTED AT THE MOMENT!!! ; exit
  fi
  let jproc=$jproc+1
done

## ======================================================================
## Restarts - OBCs 
## ======================================================================
if [ ${OBC_ACTIVE} = 'no' ] ; then 
   nomrestart=${CONFIG}_${clstep}_restart_obc.output
   expatrie_mv $nomrestart $DIR_TMP_RST $DIR_ARCHI_RST
   if [ $? -ne 0 ] ; then echo "$nomrestart non trouve" ; exit ; fi
elif [ ${OBC_ACTIVE} = 'yes' ] ; then
   nomrestart=${CONFIG}_${clstep}_restart_obc.output
   expatrie_mv $nomrestart $DIR_TMP_RST $DIR_ARCHI_RST
   if [ $? -ne 0 ] ; then echo "$nomrestart non trouve" ; exit ; fi
   
   echo 'JPP******************************************************'
   echo 'JPP - OBC_ACTIVE=yes ; OPTION NOT CODED YET!' ; #exit
   echo 'JPP - DEBUG STILL ON FOR THIS OPTION'
   echo 'JPP******************************************************'
 
fi

## ======================================================================
## Automatic re-soumission
## ======================================================================
cd $DIR_CALCU_CTL
TESTSUB=`wc $CONFIG.db`
test=` echo $TESTSUB|cut -d' ' -f1`
if [ $test -le $MAXSUB ] ; then
ord_soumet ./RAPPAT.ksh -mach pegasus -cpus 1 -mpi -cm ${MEMLIMIT} -t ${TIMELIMIT} \
                        -listing ${DIR_CALCU_CTL} -q ${PEGASUS_QUEUE}
fi


## ======================================================================
## 2D/3D outputs
## ======================================================================
## JPP MOVE OUTPUT FILES TO "STORAGE DISK"
list_dirs='DIMG LOG'
for dir in ${list_dirs} ; do
  if [ ! -d $CDIR/${CONF_NAME}/${CONFIG}-OUT/${dir} ] ; then 
     mkdir -m 755 -p $CDIR/${CONF_NAME}/${CONFIG}-OUT/${dir}
  fi
  mv -f $DIR_TMP_OUT/${dir}/* $CDIR/${CONF_NAME}/${CONFIG}-OUT/${dir}
done


## REBUILD 2D/3D OUTPUTS - JPP using basic function flio_rbld
cd $CDIR/${CONF_NAME}/${CONFIG}-OUT/DIMG

listgrids='grid_T grid_U grid_V grid_W icemod'
for grid in $listgrids ; do
  list=`ls ${CONFIG}_IN_${ndate0}_*${grid}*`
  finf=`echo $list | cut -d' ' -f1 | cut -d'_' -f5`
  outfile=${CONFIG}_IN_${ndate0}_${finf}_${grid}.nc
  echo ${DIR_ARCHI_CMB}'/'$outfile

  if [ ! -d ${DIR_ARCHI_CMB} ] ; then
    mkdir -m 755 -p ${DIR_ARCHI_CMB}
  fi


  if [ ! -s ${DIR_ARCHI_CMB}/$outfile ]; then
    #=========================================
    # begin of the rebuild script...
    #-
    # Retrieving and validation of the options
    #-
    r_v='silencious'; 
    shift $(($OPTIND-1));
    #-
    # Validate the names of the input files
    #-
    qi=0;
    for i in $list ; do
      ((qi=qi+1));
    done
    #-
    # Create the information file for the program
    #-
    let q2=$qi+1
    echo ${r_v} > tmp.$$;                  # -- Option silencious
    echo ${q2} >> tmp.$$;                  # -- # of files
    for i in $list ; do 
      echo ${i} >> tmp.$$;                 # -- Individual file names
    done
    echo ${DIR_ARCHI_CMB}/${outfile} >> tmp.$$;  # -- Output name
    #-
    # Call to the recombine routine
    #-
    cp ${DIR_CALCU_UTL}/flio_rbld .
    ./flio_rbld < tmp.$$
    #-
    # Clear
    #-
    rm -f tmp.$$ ./flio_rbld 
    # end of the rebuild script...
    #==========================================================
  else 
    echo 'JPP PP: '${DIR_ARCHI_CMB}'/'$outfile 'EXISTS'
  fi
done


## **********************************************************************
## **********************************************************************
## **********  BEGINNING OF SECTION UNUSED FOR CONCEPTS CODE   **********
## ======================================================================
## Diadct
## ======================================================================
#grep '$(prefix)key_diadct' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#if [ $? -eq 0 ] ; then
#    name=${CONFIG}_*_diadct
#    nomdct=${CONFIG}_DIADCT_${nomenclature_run}
#    mv ${DIR_TMP_DIMG}/$name ${DIR_TMP_DCT}/$nomdct
#    echo ${nomdct} $DIR_TMP_DCT $DIR_ARCHI_OUTPUT >> $DIR_TMP_OUT/list.files
#fi

## ======================================================================
## Diamoor
## ======================================================================
#moor=`grep ln_diamoor ${DIR_CALCU_CTL}/namelist-${CONFIG} | tail -1 | awk -F" " '{print $3}'`
#if [ ${moor} == ".true." ] ; then
#    cd $DIR_TMP_MOOR
#    cp $DIR_CALCU_UTL/build_nc .
#    cp $DIR_CALCU_CTL/namelistio_standard namelistio
#    sed -e "s/CONFIG/$CONFIG/"    \
#        -e "s/RDT/$rdt/"          \
#        -e "s/NDATE0/$ndate0/"    \
#        -e "s/FREQ/$nbiter/"      \
#        -e "s/FQHF/$nbiter_hf/"   \
#        -e "s/REC3D/false/"       namelistio > namelist1
#    mv namelist1 namelistio
#    ln -s $DIR_TMP_INPUT/position.moor .
#    ./build_nc 'moorings'
#    nomtar=${CONFIG}_MOORINGS_${nomenclature_run}
#    for name in *.mooring.nc ; do
#        mv ${name} ${CONFIG}_${nomenclature_run}_${name}
#    done
#    tar cvf $nomtar ${CONFIG}_${nomenclature_run}*.nc
#    echo ${nomtar} $DIR_TMP_MOOR $DIR_ARCHI_OUTPUT >> $DIR_TMP_OUT/list.files
#    rm -f *.mooring* build_nc namelistio position.moor build_nc_moorings_OK
#fi
#
## ======================================================================
## Tidal analysis results: 2d fields (elevation+barotropic velocities)
## ======================================================================
#grep '$(prefix)key_diaharm' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#if [ $? -eq 0 ] ; then
#    cd $DIR_TMP_DIMG
#    cp $DIR_CALCU_CTL/namelistio_standard namelistio
#    sed -e "s/CONFIG/$CONFIG/"    \
#        -e "s/RDT/$rdt/"          \
#        -e "s/NDATE0/$ndate0/"    \
#        -e "s/FREQ/$nbiter/"      \
#        -e "s/FQHF/$nbiter_hf/"   \
#        -e "s/REC3D/false/"       namelistio > namelist1
#    mv namelist1 namelistio
#    cp $DIR_CALCU_UTL/build_nc_harm .
#    ln -s $DIR_TMP_INPUT/coordinates.nc .
#    ./build_nc_harm
#    rm namelistio build_nc_harm coordinates.nc
#    list="T U V"
#    for var in ${list} ; do
#        name=${CONFIG}_Tidal_harmonics_grid${var}.nc
#        nametid=${CONFIG}_Tidal_harmonics_${ndate0}_${ndastpfin}_grid${var}_R${ndate0}.nc
#        if ! test -f ${name} ; then
#            echo "WARNING !!!!!!!"
#            echo ${name}
#        else
#            mv ${name} $DIR_TMP_CDF/${nametid}
#            echo ${nametid} $DIR_TMP_CDF $DIR_ARCHI_OUTPUT >> $DIR_TMP_OUT/list.files
#        fi
#    done
#fi
#
## ======================================================================
## Tidal analysis results: 3d fields (elevation+barotropic velocities)
## ======================================================================
#grep '$(prefix)key_diaharm_3d' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#if [ $? -eq 0 ] ; then
#    cd $DIR_TMP_DIMG
#    cp $DIR_CALCU_CTL/namelistio_standard namelistio
#    sed -e "s/CONFIG/$CONFIG/"    \
#        -e "s/RDT/$rdt/"          \
#        -e "s/NDATE0/$ndate0/"    \
#        -e "s/FREQ/$nbiter/"      \
#        -e "s/FQHF/$nbiter_hf/"   \
#        -e "s/REC3D/false/"       namelistio > namelist1
#    mv namelist1 namelistio
#    cp $DIR_CALCU_UTL/build_nc_harm3d .
#    ln -s $DIR_TMP_INPUT/coordinates.nc .
#    ./build_nc_harm3d 
#    rm namelistio build_nc_harm3d coordinates.nc
#    list="M2_U M2_V S2_U S2_V N2_U N2_V K1_U K1_V O1_U O1_V Q1_U Q1_V M4_U M4_V"
#    for var in ${list} ; do
#        name=${CONFIG}_${var}_Tidal_harmonics_gridT.nc
#        nametid=${CONFIG}_${var}_Tidal_harmonics_${ndate0}_${ndastpfin}_gridT_R${ndate0}.nc
#        if ! test -f ${name} ; then
#            echo "WARNING !!!!!!!"
#            echo ${name}
#        else
#            mv ${name} $DIR_TMP_CDF/${nametid}
#            echo ${nametid} $DIR_TMP_CDF $DIR_ARCHI_OUTPUT >> $DIR_TMP_OUT/list.files
#        fi
#    done
#fi
#

## ======================================================================
## Complement for 2D/3D fields 
#cd $DIR_TMP_DIMG
#list="T S U V W"
#for date in `awk -F" " '{printf "y%sm%sd%s ",substr($2,1,4),substr($2,5,2),substr($2,7,2)}' datrj.out` ; do
#    year=` echo $date | cut -c2-5 `
#    mm=` echo $date | cut -c7-8 `
#    dd=` echo $date | cut -c10-11 `
#    yymmdd="${year}${mm}${dd}"
#    yymmddp1=$(/bin/date "+%Y%m%d" --date="${yymmdd} 1 days")
#    for var in ${list} ; do
#        name=${CONFIG}_${date}_grid${var}.nc
#        nameout=${CONFIG}_1dAV_${yymmdd}_${yymmddp1}_grid${var}_R${ndate0}.nc
#        if ! test -f ${name} ; then
#            echo "WARNING !!!!!!!"
#            echo ${name} is missing
#        else
#            mv ${name} $DIR_TMP_CDF/${nameout}
#            echo ${nameout} $DIR_TMP_CDF $DIR_ARCHI_OUTPUT >> $DIR_TMP_OUT/list.files
#        fi
#    done
#done
#grep '$(prefix)key_dimgout_hf' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#if [ $? -eq 0 ] ; then
#    list="T U V"
#    for date in `awk -F" " '{printf "y%sm%sd%s ",substr($2,1,4),substr($2,5,2),substr($2,7,2)}' datrj.out` ; do
#        year=` echo $date | cut -c2-5 `
#        mm=` echo $date | cut -c7-8 `
#        dd=` echo $date | cut -c10-11 `
#        yymmdd="${year}${mm}${dd}"
#        yymmddp1=$(/bin/date "+%Y%m%d" --date="${yymmdd} 1 days")
#        for var in ${list} ; do
#            name=${CONFIG}_${date}_grid${var}HF.nc
#            nameout=${CONFIG}_1hAV_${yymmdd}_${yymmddp1}_grid${var}_R${ndate0}.nc
#            if ! test -f ${name} ; then
#                echo "WARNING !!!!!!!"
#                echo ${name} is missing
#            else
#                mv ${name} $DIR_TMP_CDF/${nameout}
#                echo ${nameout} $DIR_TMP_CDF $DIR_ARCHI_OUTPUT >> $DIR_TMP_OUT/list.files
#            fi
#        done
#    done
#fi
#
#grep '$(prefix)key_top' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#kcount=$?
#grep '$(prefix)key_pisces' $DIR_CALCU_SRC/cpp_NEATL12-T03.fcm
#((kcount=${kcount}+$?))
#if [ ${kcount} -eq 0 ] ; then
#    list="ALK BFe BSi CaCO3 DCHL DFe DIC DOC DSi Fer GOC NCHL NFe NH4 NO3 O2 PHY2 PHY PO4 POC SFe Si ZOO2 ZOO"
#    for date in `awk -F" " '{printf "y%sm%sd%s ",substr($2,1,4),substr($2,5,2),substr($2,7,2)}' datrj.out` ; do
#        year=` echo $date | cut -c2-5 `
#        mm=` echo $date | cut -c7-8 `
#        dd=` echo $date | cut -c10-11 `
#        yymmdd="${year}${mm}${dd}"
#        yymmddp1=$(/bin/date "+%Y%m%d" --date="${yymmdd} 1 days")
#        for var in ${list} ; do
#            name=${CONFIG}_${date}_grid${var}.nc
#            nameout=${CONFIG}_1dAV_${yymmdd}_${yymmddp1}_grid${var}_R${ndate0}.nc
#            if ! test -f ${name} ; then
#                echo "WARNING !!!!!!!"
#                echo ${name} is missing
#            else
#                mv ${name} $DIR_TMP_CDF/${nameout}
#                echo ${nameout} $DIR_TMP_CDF $DIR_ARCHI_OUTPUT >> $DIR_TMP_OUT/list.files
#            fi
#        done
#    done
#fi
#
#rm $DIR_TMP_DIMG/datrj.out
#
## **********      END OF SECTION UNUSED FOR CONCEPTS CODE     **********
## **********************************************************************
## **********************************************************************




## ======================================================================
## Submission of the ARCHIV script
## ======================================================================
if [ $archive = "true" ] ; then
    cd $DIR_CALCU_CTL
    qsub ARCHIV.ksh
else
    rm $DIR_TMP_OUT/list.files
fi
