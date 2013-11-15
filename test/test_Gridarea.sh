#!/bin/sh
#
CDODEBUG=0
#
if [ "$CDODEBUG" = 0 ]; then CDO="$CDO -s"; fi
CDOOUT=cout
CDOERR=cerr
REFVAL="12.5663706"
GRIDS="global_5 global_2 global_1 global_.5 n32 n80 n160"
RSTAT=0;
PLANET_RADIUS=1
export PLANET_RADIUS
#
for GRID in $GRIDS; do
  GLOBAREA=`$CDO outputf,%10.7f -fldsum -gridarea -random,$GRID`
#  echo "$GRID: >$GLOBAREA< >$REFVAL<"
  if [ "$GLOBAREA" != "$REFVAL" ]; then RSTAT=`expr $RSTAT + 1`; fi
done
#
rm -f $CDOOUT $CDOERR
unset PLANET_RADIUS
#
if [ "$CDODEBUG" = 1 ]; then
  echo "rstat: $RSTAT"
fi
#
exit $RSTAT
