#!/bin/sh
#
#CDO=cdo
#DATAPATH=data
#
CDODEBUG=0
#
if [ "$CDODEBUG" = 0 ]; then CDO="$CDO -s"; fi
CDOOUT=cout
CDOERR=cerr
FORMAT="-f srv -b 32"
GRIDS="n16 n32"
RSTAT=0;
#
IFILE=$DATAPATH/bathy4.grb
#
for GRIDTYPE in " " "-setgridtype,curvilinear" "-setgridtype,unstructured"; do
  for GRID in $GRIDS; do
# remaplaf: sort could give different results"
    RMODS="bil bic nn con"
    if [ "$GRIDTYPE" = "-setgridtype,unstructured" ]; then RMODS="nn con"; fi
    for RMOD in $RMODS; do
      OFILE=${GRID}_${RMOD}
      RFILE=$DATAPATH/${OFILE}_ref
      $CDO $FORMAT remap${RMOD},$GRID $GRIDTYPE $IFILE ${OFILE} > $CDOOUT 2> $CDOERR
      if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
      if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
      $CDO diff $OFILE $RFILE > $CDOOUT 2> $CDOERR
      if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
      if [ -s $CDOOUT ] ; then RSTAT=`expr $RSTAT + 1`; fi
      if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
      rm -f $OFILE
    done
  done
done
#
rm -f $CDOOUT $CDOERR
#
if [ "$CDODEBUG" = 1 ]; then
  echo "rstat: $RSTAT"
fi
#
exit $RSTAT
