#!/bin/sh
#
CDODEBUG=0
#
if [ "$CDODEBUG" = 0 ]; then CDO="$CDO -s"; fi
CDOOUT=cout
CDOERR=cerr
FORMAT=""
RSTAT=0;
#
IFILE=$DATAPATH/t21_geosp_tsurf.grb
RFILE=$DATAPATH/gp2sp_ref
OFILE=gp2sp_res
$CDO $FORMAT gp2sp $IFILE $OFILE > $CDOOUT 2> $CDOERR
if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
$CDO diff $OFILE $RFILE > $CDOOUT 2> $CDOERR
if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
if [ -s $CDOOUT ] ; then RSTAT=`expr $RSTAT + 1`; fi
if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
rm -f $OFILE
#
IFILE=$DATAPATH/gp2sp_ref
RFILE=$DATAPATH/sp2gp_ref
OFILE=sp2gp_res
$CDO $FORMAT sp2gp $IFILE $OFILE > $CDOOUT 2> $CDOERR
if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
$CDO diff $OFILE $RFILE > $CDOOUT 2> $CDOERR
if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
if [ -s $CDOOUT ] ; then RSTAT=`expr $RSTAT + 1`; fi
if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
rm -f $OFILE
#
IFILE=$DATAPATH/t21_geosp_tsurf.grb
RFILE=$DATAPATH/gp2spl_ref
OFILE=gp2spl_res
$CDO $FORMAT gp2spl $IFILE $OFILE > $CDOOUT 2> $CDOERR
if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
$CDO diff $OFILE $RFILE > $CDOOUT 2> $CDOERR
if [ $? != 0 ]    ; then RSTAT=`expr $RSTAT + 1`; fi
if [ -s $CDOOUT ] ; then RSTAT=`expr $RSTAT + 1`; fi
if [ "$CDODEBUG" = 1 ]; then cat $CDOOUT $CDOERR; fi
rm -f $OFILE
#
IFILE=$DATAPATH/gp2spl_ref
RFILE=$DATAPATH/sp2gpl_ref
OFILE=sp2gpl_res
$CDO $FORMAT sp2gpl $IFILE $OFILE > $CDOOUT 2> $CDOERR
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
