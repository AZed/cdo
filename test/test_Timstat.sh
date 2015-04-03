#!/bin/sh
#
CDODEBUG=0
#
if [ "$CDODEBUG" = 0 ]; then CDO="$CDO -s"; fi
CDOOUT=cout
CDOERR=cerr
STATS="min max sum avg mean std std1 var var1"
RSTAT=0;
#
IFILE=$DATAPATH/ts_mm_5years
#
for STAT in $STATS; do
  RFILE=$DATAPATH/tim${STAT}_ref
  OFILE=tim${STAT}_res
  $CDO tim${STAT} $IFILE $OFILE > $CDOOUT 2> $CDOERR
  if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
  if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
  $CDO diff $OFILE $RFILE > $CDOOUT 2> $CDOERR
  if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
  if [ -s $CDOOUT ] ; then RSTAT=`expr $RSTAT + 1`; fi
  if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
  rm -f $OFILE
done
#
rm -f $CDOOUT $CDOERR
#
if [ "$CDODEBUG" = 1 ]; then
  echo "rstat: $RSTAT"
fi
#
exit $RSTAT
