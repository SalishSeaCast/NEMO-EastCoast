#!/bin/ksh
## ======================================================================
##            SCRIPT TO SUCCESSIVELY SUBMIT JOBS ON BALTIC
##            $CONF.ksh
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

# Fred...
typeset -Z5 irstrt
typeset -Z5 bustime
typeset -Z5 ncpus

set -x
echo `pwd`
. ./environment.ksh

## ======================================================================
## Execution of OPA
## ======================================================================
echo '*******************'
echo '* RUN * RUN * RUN *'
echo '*******************'

cd $DIR_TMP_RUN


# Prepare dirs for output files for each procs
ncpus=00001
while [ ${ncpus} -le $NBPROC ] ; do
  mkdir ${ncpus}
  ncpus=`expr ${ncpus} + 1`
done

ulimit -s unlimited  # increase memory
ulimit -t ${TIMELIMIT}
rumpirun.openmpi -np ${NBPROC} ./opa


## ======================================================================
## We save all the log files and we test if the run is OK or not
## ======================================================================
## Move log files 
mv ocean.output* ${DIR_TMP_LOG}/.
mv time.step ${DIR_TMP_LOG}/.
mv layout.dat ${DIR_TMP_LOG}/.
mv ${CONFIG}.db ${DIR_TMP_LOG}/.
mv namelist* ${DIR_TMP_LOG}/.
mv solver.stat ${DIR_TMP_LOG}/.

## Move Restarts and outputs out of the execution directory
## (Archives and definitive move done in POST_TREATMENT.ksh)
mv ${CONFIG}*_restart_obc.output ${DIR_TMP_RST}/.  #- Restarts OBCs
mv *_restart_*.nc ${DIR_TMP_RST}/.		   #- Restarts ocean and sea ice 
mv  ${CONFIG}*.nc ${DIR_TMP_DIMG}/.                #- Model output
if [ -f  output.init_0001.nc ] ; then		   #- Initial conditions if active
  mv output.init_00*.nc ${DIR_TMP_OUT}
fi

## Clean DIRRUN...
#rm -fr ATMDATA 
# Remove dirs for output files for each procs
ncpus=00001
while [ ${ncpus} -le $NBPROC ] ; do
  rm -fr ${ncpus}
  ncpus=`expr ${ncpus} + 1`
done


echo CALL TO POST-TREATMENT.ksh  !!!!!!!!!!!
echo ~~~~~~~~~~~~~~~~~~~~~
cd $DIR_CALCU_CTL
ord_soumet ./POST-TREATMENT.ksh -mach pegasus -cpus 1 -mpi -cm ${MEMLIMIT} -t ${TIMELIMIT} \
                                -listing ${DIR_CALCU_CTL} -q ${PEGASUS_QUEUE}

## ======================================================================
## End of the $SCONF script
## ======================================================================
