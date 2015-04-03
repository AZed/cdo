#include <stdio.h>		/* for NULL */
#include <math.h>		/* for modf(), floor(), log10(), ceil() */
#include <float.h>		/* for DBL_MANT_DIG */

#include "cdi.h"  		/* CALENDAR_ */

#undef 	 ABS
#define	 ABS(a)		((a) < 0 ? -(a) : (a))


static int month_360[12] = {30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30};
static int month_365[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static int month_366[12] = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};



void decode_date(int date, int *year, int *month, int *day);
int encode_date(int year, int month, int day);

void decode_time(int time, int *hour, int *minute);
int encode_time(int hour, int minute);


/*
 * The following two functions convert between Julian day number and
 * Gregorian/Julian dates (Julian dates are used prior to October 15,
 * 1582; Gregorian dates are used after that).  Julian day number 0 is
 * midday, January 1, 4713 BCE.  The Gregorian calendar was adopted
 * midday, October 15, 1582.
 *
 * Author: Robert Iles, March 1994
 *
 * C Porter: Steve Emmerson, October 1995
 *
 * Original: http://www.nag.co.uk:70/nagware/Examples/calendar.f90
 *
 * There is no warranty on this code.
 */


/*
 * Convert a Julian day number to a Gregorian/Julian date.
 */
static void julday_to_gregdate(
    unsigned long	julday,		/* Julian day number to convert */
    int			*year,		/* Gregorian year (out) */
    int			*month,		/* Gregorian month (1-12) (out) */
    int			*day		/* Gregorian day (1-31) (out) */
    )
{
#if INT_MAX <= 0X7FFF
    long	ja, jb, jd;
#else
    int		ja, jb, jd;
#endif
    int		jc;
    int		je, iday, imonth, iyear;
    double	xc;

    if (julday < 2299161)
	ja = julday;
    else
    {
      int ia = (int) (((julday - 1867216) - 0.25) / 36524.25);

	ja = julday + 1 + ia - (int)(0.25 * ia);
    }

    jb = ja + 1524;
    xc = ((jb - 2439870) - 122.1) / 365.25;
    jc = (int) (6680.0 + xc);
    jd = 365 * jc + (int)(0.25 * jc);
    je = (int)((jb - jd) / 30.6001);

    iday = (int)(jb - jd - (int)(30.6001 * je));

    imonth = je - 1;
    if (imonth > 12)
	imonth -= 12;

    iyear = jc - 4715;
    if (imonth > 2)
	iyear -= 1;
    if (iyear <= 0)
	iyear -= 1;

    *year = iyear;
    *month = imonth;
    *day = iday;
}


/*
 * Convert a Gregorian/Julian date to a Julian day number.
 *
 * The Gregorian calendar was adopted midday, October 15, 1582.
 */
static unsigned long gregdate_to_julday(
    int		year,	/* Gregorian year */
    int		month,	/* Gregorian month (1-12) */
    int		day	/* Gregorian day (1-31) */
    )
{
#if INT_MAX <= 0X7FFF
    long		igreg = 15 + 31 * (10 + (12 * 1582));
    long		iy;	/* signed, origin 0 year */
    long		ja;	/* Julian century */
    long		jm;	/* Julian month */
    long		jy;	/* Julian year */
#else
    int			igreg = 15 + 31 * (10 + (12 * 1582));
    int			iy;	/* signed, origin 0 year */
    int			ja;	/* Julian century */
    int			jm;	/* Julian month */
    int			jy;	/* Julian year */
#endif
    unsigned long	julday;	/* returned Julian day number */

    /*
     * Because there is no 0 BC or 0 AD, assume the user wants the start of 
     * the common era if they specify year 0.
     */
    if (year == 0)
	year = 1;

    iy = year;
    if (year < 0)
	iy++;
    if (month > 2)
    {
	jy = iy;
	jm = month + 1;
    }
    else
    {
	jy = iy - 1;
	jm = month + 13;
    }

    /*
     *  Note: SLIGHTLY STRANGE CONSTRUCTIONS REQUIRED TO AVOID PROBLEMS WITH
     *        OPTIMISATION OR GENERAL ERRORS UNDER VMS!
     */
    julday = day + (int)(30.6001 * jm);
    if (jy >= 0)
    {
	julday += 365 * jy;
	julday += (unsigned long) (0.25 * jy);
    }
    else
    {
	double		xi = 365.25 * jy;

	if ((int)xi != xi)
	    xi -= 1;
	julday += (int)xi;
    }
    julday += 1720995;

    if (day + (31* (month + (12 * iy))) >= igreg)
    {
	ja = jy/100;
	julday -= ja;
	julday += 2;
	julday += ja/4;
    }

    return julday;
}


/*
 * Encode a date as a double-precision value.
 */
static double utencdate(int year, int month, int day)
{
    return ((long)gregdate_to_julday(year, month, day) - 
	    (long)gregdate_to_julday(2001, 1, 1)) * 86400.0;
}


/*
 * Encode a time as a double-precision value.
 */
static double utencclock(int hours, int minutes, double seconds)
{
    return (hours*60 + minutes)*60 + seconds;
}


/*
 * Decompose a value into a set of values accounting for uncertainty.
 */
static void decomp(
       double	value,
       double	uncer,		/* >= 0 */
       int	nbasis,
       double	*basis,		/* all values > 0 */
       double	*count
    )
{
    int		i;

    for (i = 0; i < nbasis; i++)
    {
	double	r = fmod(value, basis[i]);	/* remainder */

	/* Adjust remainder to minimum magnitude. */
	if (ABS(2*r) > basis[i])
	    r += r > 0
		    ? -basis[i]
		    :  basis[i];

	if (ABS(r) <= uncer)
	{
	    /* The value equals a basis multiple within the uncertainty. */
	    double	half = value < 0 ? -basis[i]/2 : basis[i]/2;
	    modf((value+half)/basis[i], count+i);
	    break;
	}

	/* value = basis[i] * modf((value+0.0000001)/basis[i], count+i); */
	value = basis[i] * modf(value/basis[i], count+i);
    }

    if (i >= nbasis) {
	count[--i] += value;
    }
    else {
	for (i++; i < nbasis; i++)
	    count[i] = 0;
    }
}


/*
 * Decode a time from a double-precision value.
 */
static void dectime(
	double	value,
	int	*year,
	int	*month,
	int	*day,
	int	*hour,
	int	*minute,
	float	*second)
{
    long	days;
    long	hours;
    long	minutes;
    double	seconds;
    double	uncer;		/* uncertainty of input value */
    typedef union
    {
	double	    vec[7];
	struct
	{
	    double	days;
	    double	hours12;
	    double	hours;
	    double	minutes10;
	    double	minutes;
	    double	seconds10;
	    double	seconds;
	}	    ind;
    } Basis;
    Basis	counts;
    static Basis	basis = {{86400, 43200, 3600, 600, 60, 10, 1}};

    uncer = ldexp(value < 0 ? -value : value, -DBL_MANT_DIG);

    days = (long) floor(value/basis.ind.days);
    value -= days * basis.ind.days;		/* make positive excess */

    decomp(value, uncer, sizeof(basis.vec)/sizeof(basis.vec[0]),
	   basis.vec, counts.vec);

    days += (long) counts.ind.days;
    hours = (int)counts.ind.hours12 * 12 + (int)counts.ind.hours;
    minutes = (int)counts.ind.minutes10 * 10 + (int)counts.ind.minutes;
    seconds = (int)counts.ind.seconds10 * 10 + counts.ind.seconds;

    seconds = (float)seconds;

    if ((float)seconds >= 60) {
	seconds -= 60;
	if (++minutes >= 60) {
	    minutes -= 60;
	    if (++hours >= 24) {
		hours -= 24;
		days++;
	    }
	}
    }

    *second = seconds;
    *minute = minutes;
    *hour = hours;

    julday_to_gregdate(gregdate_to_julday(2001, 1, 1) + days, year, month, day);
}

/*
 * Convert a temporal value into XXX days/year date and time.
 *
 */
static void dectimeXXX(int dpy, double value, int *year, int *month,
		       int *day, int *hour, int *minute, float *second)
{
    int i = 0;
    int *dpm = NULL;
    long	days;
    int  	hours;
    int 	minutes;
    double	seconds;
    double	uncer;		/* uncertainty of input value */
    typedef union
    {
	double	    vec[7];
	struct
	{
	    double	days;
	    double	hours12;
	    double	hours;
	    double	minutes10;
	    double	minutes;
	    double	seconds10;
	    double	seconds;
	}	    ind;
    } Basis;
    Basis	counts;
    static Basis	basis;

    basis.ind.days      = 86400;
    basis.ind.hours12   = 43200;
    basis.ind.hours     = 3600;
    basis.ind.minutes10 = 600;
    basis.ind.minutes   = 60;
    basis.ind.seconds10 = 10;
    basis.ind.seconds   = 1;

    uncer = ldexp(value < 0 ? -value : value, -DBL_MANT_DIG);

    days = (long) floor(value/86400.0);
    value -= days * 86400.0;		/* make positive excess */

    decomp(value, uncer, sizeof(basis.vec)/sizeof(basis.vec[0]),
	   basis.vec, counts.vec);

    days   += (long) counts.ind.days;
    hours   = (int)counts.ind.hours12   * 12 + (int)counts.ind.hours;
    minutes = (int)counts.ind.minutes10 * 10 + (int)counts.ind.minutes;
    seconds = (int)counts.ind.seconds10 * 10 + counts.ind.seconds;

    *second = seconds;
    *minute = minutes;
    *hour   = hours;

    *year = (days-1) / dpy;

    days -= (*year*dpy);

    if      ( dpy == 360 ) dpm = month_360;
    else if ( dpy == 365 ) dpm = month_365;
    else if ( dpy == 366 ) dpm = month_366;

    if ( dpm )
      for ( i = 0; i < 12; i++ )
	{
	    if ( days > dpm[i] ) days -= dpm[i];
	    else break;
	}

    *month = i + 1;

    *day = days;
}


/*
 * Encode a date as a double-precision value.
 */
static double encdateXXX(int dpy, int year, int month, int day)
{
  int i;
  int *dpm = NULL;
  double rval;

  rval = dpy * year + day;

  if      ( dpy == 360 ) dpm = month_360;
  else if ( dpy == 365 ) dpm = month_365;
  else if ( dpy == 366 ) dpm = month_366;
  
  if ( dpm ) for ( i = 0; i < month-1; i++ ) rval += dpm[i];

  return (rval * 86400.0);
}


/*
 * Convert a Gregorian/Julian date and time into a temporal value.
 *
 */
double encode_julval(int dpy, int date, int time)
{
  int year, month, day, hour, minute;
  double second = 0;
  double julval;

  decode_date(date, &year, &month, &day);
  decode_time(time, &hour, &minute);
  
  if ( dpy == 360 || dpy == 365 || dpy == 366 )
    julval = encdateXXX(dpy, year, month, day) + utencclock(hour, minute, second);
  else
    julval = utencdate(year, month, day) + utencclock(hour, minute, second);

  return (julval);
}


/*
 * Convert a temporal value into UTC Gregorian/Julian date and time.
 *
 */
void decode_julval(int dpy, double value, int *date, int *time)
{
  int year, month, day, hour, minute;
  float	sec;

  if ( dpy == 360 || dpy == 365 || dpy == 366 )
    dectimeXXX(dpy, value, &year, &month, &day, &hour, &minute, &sec);
  else
    dectime(value, &year, &month, &day, &hour, &minute, &sec);
  /*
  fprintf(stdout, "%d %d %d %d %d %g\n", year, month, day, hour, minute, (double)sec);
  */
  *date = encode_date(year, month, day);
  *time = encode_time(hour, minute);
}
