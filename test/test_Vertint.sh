#!/bin/sh
#
CDODEBUG=0
#
if [ "$CDODEBUG" = 0 ]; then CDO="$CDO -s"; fi
CDOOUT=cout
CDOERR=cerr
FORMAT="-f srv -b 32"
RSTAT=0;
#
IFILE=$DATAPATH/hl_l19.grb
#
RFILE=$DATAPATH/ml2pl_ref
OFILE=ml2pl_res
$CDO $FORMAT ml2pl,92500,85000,50000,20000 $IFILE $OFILE > $CDOOUT 2> $CDOERR
if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
$CDO diff $OFILE $RFILE > $CDOOUT 2> $CDOERR
if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
if [ -s $CDOOUT ] ; then RSTAT=`expr $RSTAT + 1`; fi
if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
rm -f $OFILE
#
rm -f $CDOOUT $CDOERR
#
if [ "$CDODEBUG" = 1 ]; then
  echo "rstat: $RSTAT"
fi
#
exit $RSTAT
