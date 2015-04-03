#!/bin/sh
#
#CDO=cdo
#DATAPATH=data
#
CDODEBUG=0
#
if [ "$CDO_TEST_DEBUG" = 1 ]; then CDODEBUG=1; fi
#
if [ "$CDODEBUG" = 0 ]; then CDO="$CDO -s"; fi
CDOOUT=cout
CDOERR=cerr
FORMAT="-f srv -b 32"
RSTAT=0;
#
IFILE=$DATAPATH/t21_geosp_tsurf.grb
#
RFILE=catdata_ref
OFILE=catdata
#
cp $IFILE ${RFILE}
chmod u+w ${RFILE}
cat $IFILE >> ${RFILE}
#
rm -f ${OFILE}
$CDO cat $IFILE ${OFILE}
$CDO cat $IFILE ${OFILE}
#
$CDO diff $OFILE $RFILE > $CDOOUT 2> $CDOERR
if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
if [ -s $CDOOUT ] ; then RSTAT=`expr $RSTAT + 1`; fi
if [ -s $CDOERR ] ; then RSTAT=`expr $RSTAT + 1`; fi
if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
rm -f $OFILE $RFILE
#
rm -f $CDOOUT $CDOERR
#
if [ "$CDODEBUG" = 1 ]; then
  echo "rstat: $RSTAT"
fi
#
exit $RSTAT
