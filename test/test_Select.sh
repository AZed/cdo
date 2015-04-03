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
IFILE=$DATAPATH/pl_data.grb
#
TNUM=0
#
for SELECT in "code=130,152" "code=130,152,level=9000,90" "code=130,152,level=9000,90,timestep=2,3,5" "code=130" "level=90000"; do
  TNUM=`expr $TNUM + 1`

  RFILE=$DATAPATH/select${TNUM}_ref
  OFILE=select${TNUM}_res
  $CDO select,${SELECT} $IFILE $OFILE > $CDOOUT 2> $CDOERR
  if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
  if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
  $CDO diff $OFILE $RFILE > $CDOOUT 2> $CDOERR
  if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
  if [ -s $CDOOUT ] ; then RSTAT=`expr $RSTAT + 1`; fi
  if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
  rm -f $OFILE
#
  rm -f $CDOOUT $CDOERR
done
#
TNUM=0
#
for DELETE in "code=129" "code=129,130,level=90000,900" "code=129,130,level=90000,900,timestep=1,4"; do
  TNUM=`expr $TNUM + 1`

  RFILE=$DATAPATH/select${TNUM}_ref
  OFILE=delete${TNUM}_res
  $CDO delete,${DELETE} $IFILE $OFILE > $CDOOUT 2> $CDOERR
  if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
  if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
  $CDO diff $OFILE $RFILE > $CDOOUT 2> $CDOERR
  if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
  if [ -s $CDOOUT ] ; then RSTAT=`expr $RSTAT + 1`; fi
  if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
  rm -f $OFILE
#
  rm -f $CDOOUT $CDOERR
done
#
if [ "$CDODEBUG" = 1 ]; then
  echo "rstat: $RSTAT"
fi
#
exit $RSTAT
