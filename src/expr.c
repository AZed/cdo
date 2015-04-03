#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "field.h"
#include "expr.h"
#include "expr_yacc.h"


static double f_abs(double x)  { return (fabs(x));  }
static double f_int(double x)  { return ((int)(x)); }
static double f_nint(double x) { return (NINT(x));  }
static double f_sqr(double x)  { return (x*x);      }

typedef struct {
  int type;
  char *name;                      /* function name            */
  double (*func)(double);          /* pointer to function      */
}
func_t;

static func_t fun_sym_tbl[] =
{
  /* scalar functions */
  {0, "abs",   f_abs},
  {0, "int",   f_int},
  {0, "nint",  f_nint},
  {0, "sqr",   f_sqr},
  {0, "sqrt",  sqrt},
  {0, "exp",   exp},
  {0, "log",   log},
  {0, "log10", log10},
  {0, "sin",   sin},
  {0, "cos",   cos},
  {0, "tan",   tan},
  {0, "asin",  asin},
  {0, "acos",  acos},
  {0, "atan",  atan},

  /* array functions
  {1, "min",   min},
  {1, "max",   max},
  {1, "sum",   sum},
  {1, "avg",   avg},
  {1, "mean",  mean},
  {1, "std",   std},
  {1, "var",   var},
  */
};

static int NumFunc = sizeof(fun_sym_tbl) / sizeof(fun_sym_tbl[0]);

static
nodeType *expr_con_con(int oper, nodeType *p1, nodeType *p2)
{
  nodeType *p;

  p = (nodeType *) malloc(sizeof(nodeType));

  p->type = typeCon;

  switch ( oper )
    {
    case '+':  p->u.con.value = p1->u.con.value + p2->u.con.value; break;
    case '-':  p->u.con.value = p1->u.con.value - p2->u.con.value; break;
    case '*':  p->u.con.value = p1->u.con.value * p2->u.con.value; break;
    case '/':  p->u.con.value = p1->u.con.value / p2->u.con.value; break;
    case '^':  p->u.con.value = pow(p1->u.con.value, p2->u.con.value); break;
    default:   cdoAbort("%s: operator %c unsupported!", __func__, oper);
    }

  return (p);
}

static
nodeType *expr_con_var(int oper, nodeType *p1, nodeType *p2)
{
  nodeType *p;
  long ngp, i;
  long nlev;
  int nmiss;
  int gridID, zaxisID;
  double missval1, missval2;

  gridID   = p2->gridID;
  zaxisID  = p2->zaxisID;
  nmiss    = p2->nmiss;
  missval1 = p2->missval;
  missval2 = p2->missval;

  ngp  = gridInqSize(gridID);
  nlev = zaxisInqSize(zaxisID);

  p = (nodeType *) malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->tmpvar   = 1;
  p->u.var.nm = strdupx("tmp");
  p->gridID   = gridID;
  p->zaxisID  = zaxisID;
  p->missval  = missval1;

  p->data = (double *) malloc(ngp*nlev*sizeof(double));

  switch ( oper )
    {
    case '+':
      if ( nmiss > 0 )
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = ADD(p1->u.con.value, p2->data[i]);
	}
      else
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = p1->u.con.value + p2->data[i];
	}
      break;
    case '-':
      if ( nmiss > 0 )
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = SUB(p1->u.con.value, p2->data[i]);
	}
      else
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = p1->u.con.value - p2->data[i];
	}
      break;
    case '*':
      if ( nmiss > 0 )
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = MUL(p1->u.con.value, p2->data[i]);
	}
      else
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = p1->u.con.value * p2->data[i];
	}
      break;
    case '/':
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = DIV(p1->u.con.value, p2->data[i]);
	}
      break;
    case '^':
      if ( nmiss > 0 )
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = POW(p1->u.con.value, p2->data[i]);
	}
      else
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = pow(p1->u.con.value, p2->data[i]);
	}
      break;
    default:
      cdoAbort("%s: operator %c unsupported!", __func__, oper);
    }

  nmiss = 0;
  for ( i = 0; i < ngp*nlev; i++ )
    if ( DBL_IS_EQUAL(p->data[i], missval1) ) nmiss++;

  p->nmiss = nmiss;

  if ( p2->tmpvar ) free(p2->data);

  return (p);
}

static
nodeType *expr_var_con(int oper, nodeType *p1, nodeType *p2)
{
  nodeType *p;
  long ngp, i;
  long nlev;
  int nmiss;
  int gridID, zaxisID;
  double missval1, missval2;

  gridID   = p1->gridID;
  zaxisID  = p1->zaxisID;
  nmiss    = p1->nmiss;
  missval1 = p1->missval;
  missval2 = p1->missval;

  ngp  = gridInqSize(gridID);
  nlev = zaxisInqSize(zaxisID);

  p = (nodeType *) malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->tmpvar   = 1;
  p->u.var.nm = strdupx("tmp");
  p->gridID   = gridID;
  p->zaxisID  = zaxisID;
  p->missval  = missval1;

  p->data = (double *) malloc(ngp*nlev*sizeof(double));

  switch ( oper )
    {
    case '+':
      if ( nmiss > 0 )
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = ADD(p1->data[i], p2->u.con.value);
	}
      else
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = p1->data[i] + p2->u.con.value;
	}
      break;
    case '-':
      if ( nmiss > 0 )
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = SUB(p1->data[i], p2->u.con.value);
	}
      else
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = p1->data[i] - p2->u.con.value;
	}
      break;
    case '*':
      if ( nmiss > 0 )
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = MUL(p1->data[i], p2->u.con.value);
	}
      else
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = p1->data[i] * p2->u.con.value;
	}
      break;
    case '/':
      if ( nmiss > 0 || IS_EQUAL(p2->u.con.value, 0) )
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = DIV(p1->data[i], p2->u.con.value);
	}
      else
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = p1->data[i] / p2->u.con.value;
	}
      break;
    case '^':
      if ( nmiss > 0 )
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = POW(p1->data[i], p2->u.con.value);
	}
      else
	{
	  for ( i = 0; i < ngp*nlev; i++ )
	    p->data[i] = pow(p1->data[i], p2->u.con.value);
	}
      break;
    default:
      cdoAbort("%s: operator %c unsupported!", __func__, oper);
    }

  nmiss = 0;
  for ( i = 0; i < ngp*nlev; i++ )
    if ( DBL_IS_EQUAL(p->data[i], missval1) ) nmiss++;

  p->nmiss = nmiss;

  if ( p1->tmpvar ) free(p1->data);

  return (p);
}

static
nodeType *expr_var_var(int oper, nodeType *p1, nodeType *p2)
{
  nodeType *p;
  long ngp, ngp1, ngp2, i;
  long nlev, nlev1, nlev2, k;
  long loff1, loff2;
  int nmiss, nmiss1, nmiss2;
  double missval1, missval2;

  nmiss1   = p1->nmiss;
  nmiss2   = p2->nmiss;
  missval1 = p1->missval;
  missval2 = p2->missval;

  ngp1 = gridInqSize(p1->gridID);
  ngp2 = gridInqSize(p2->gridID);

  if ( ngp1 != ngp2 )
    cdoAbort("number of grid points differ. ngp1 = %d, ngp2 = %d", ngp1, ngp2);

  ngp = ngp1;

  nlev1 = zaxisInqSize(p1->zaxisID);
  nlev2 = zaxisInqSize(p2->zaxisID);

  p = (nodeType *) malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->tmpvar   = 1;
  p->u.var.nm = strdupx("tmp");

  if ( nlev1 > nlev2 )
    {
      nlev = nlev1;
      p->gridID  = p1->gridID;
      p->zaxisID = p1->zaxisID;
      p->missval = p1->missval;
      if ( nlev2 != 1 ) cdoAbort("nlev2 = %d must be 1!", nlev2);
    }
  else if ( nlev2 > nlev1 )
    {
      nlev = nlev2;
      p->gridID  = p2->gridID;
      p->zaxisID = p2->zaxisID;
      p->missval = p2->missval;
      if ( nlev1 != 1 ) cdoAbort("nlev1 = %d must be 1!", nlev1);
    }
  else
    {
      nlev = nlev1;
      p->gridID  = p1->gridID;
      p->zaxisID = p1->zaxisID;
      p->missval = p1->missval;
    }

  p->data = (double *) malloc(ngp*nlev*sizeof(double));

  for ( k = 0; k < nlev; k++ )
    {
      if ( nlev1 == 1 ) loff1 = 0;
      else              loff1 = k*ngp;

      if ( nlev2 == 1 ) loff2 = 0;
      else              loff2 = k*ngp;

      switch ( oper )
	{
	case '+':
	  if ( nmiss1 > 0 || nmiss2 > 0 )
	    {
	      for ( i = 0; i < ngp; i++ )
		p->data[i+k*ngp] = ADD(p1->data[i+loff1], p2->data[i+loff2]);
	    }
	  else
	    {
	      for ( i = 0; i < ngp; i++ )
		p->data[i+k*ngp] = p1->data[i+loff1] + p2->data[i+loff2];
	    }
	  break;
	case '-':
	  if ( nmiss1 > 0 || nmiss2 > 0 )
	    {
	      for ( i = 0; i < ngp; i++ )
		p->data[i+k*ngp] = SUB(p1->data[i+loff1], p2->data[i+loff2]);
	    }
	  else
	    {
	      for ( i = 0; i < ngp; i++ )
		p->data[i+k*ngp] = p1->data[i+loff1] - p2->data[i+loff2];
	    }
	  break;
	case '*':
	  if ( nmiss1 > 0 || nmiss2 > 0 )
	    {
	      for ( i = 0; i < ngp; i++ )
		p->data[i+k*ngp] = MUL(p1->data[i+loff1], p2->data[i+loff2]);
	    }
	  else
	    {
	      for ( i = 0; i < ngp; i++ )
		p->data[i+k*ngp] = p1->data[i+loff1] * p2->data[i+loff2];
	    }
	  break;
	case '/':
	  if ( nmiss1 > 0 || nmiss2 > 0 )
	    {
	      for ( i = 0; i < ngp; i++ )
		p->data[i+k*ngp] = DIV(p1->data[i+loff1], p2->data[i+loff2]);
	    }
	  else
	    {
	      for ( i = 0; i < ngp; i++ )
		p->data[i+k*ngp] = p1->data[i+loff1] / p2->data[i+loff2];
	    }
	  break;
	case '^':
	  if ( nmiss1 > 0 || nmiss2 > 0 )
	    {
	      for ( i = 0; i < ngp; i++ )
		p->data[i+k*ngp] = POW(p1->data[i+loff1], p2->data[i+loff2]);
	    }
	  else
	    {
	      for ( i = 0; i < ngp; i++ )
		p->data[i+k*ngp] = pow(p1->data[i+loff1], p2->data[i+loff2]);
	    }
	  break;
	default:
	  cdoAbort("%s: operator %c unsupported!", __func__, oper);
	}
    }

  nmiss = 0;
  for ( i = 0; i < ngp*nlev; i++ )
    if ( DBL_IS_EQUAL(p->data[i], missval1) ) nmiss++;

  p->nmiss = nmiss;

  if ( p1->tmpvar ) free(p1->data);
  if ( p2->tmpvar ) free(p2->data);

  return (p);
}

static
void ex_copy(nodeType *p2, nodeType *p1)
{
  long ngp, ngp1, ngp2, i;
  long nlev;

  if ( cdoVerbose )
    printf("\tcopy %s\n", p1->u.var.nm);

  ngp1 = gridInqSize(p1->gridID);
  ngp2 = gridInqSize(p2->gridID);

  ngp = ngp2;
  nlev = zaxisInqSize(p2->zaxisID);

  for ( i = 0; i < ngp*nlev; i++ )
    p2->data[i] = p1->data[i];

  p2->missval = p1->missval;
  p2->nmiss   = p1->nmiss;
}

static
nodeType *expr(int oper, nodeType *p1, nodeType *p2)
{
  nodeType *p = NULL;

  if ( p1->type == typeVar && p2->type == typeVar )
    {
      p = expr_var_var(oper, p1, p2);
      if ( cdoVerbose )
	printf("\t%s %c %s\n", p1->u.var.nm, oper, p2->u.var.nm);
    }
  else if ( p1->type == typeCon && p2->type == typeCon )
    {
      p = expr_con_con(oper, p1, p2);
      if ( cdoVerbose )
	printf("\t%g %c %g\n", p1->u.con.value, oper, p2->u.con.value);
    }
  else if ( p1->type == typeVar && p2->type == typeCon )
    {
      p = expr_var_con(oper, p1, p2);
      if ( cdoVerbose )
	printf("\t%s %c %g\n", p1->u.var.nm, oper, p2->u.con.value);
    }
  else if ( p1->type == typeCon && p2->type == typeVar )
    {
      p = expr_con_var(oper, p1, p2);
      if ( cdoVerbose )
	printf("\t%g %c %s\n", p1->u.con.value, oper, p2->u.var.nm);
    }
  else
    cdoAbort("Internal problem!");

  return (p);
}

static
nodeType *ex_fun_con(char *fun, nodeType *p1)
{
  nodeType *p;
  int i;
  int funcID = -1;

  p = (nodeType *) malloc(sizeof(nodeType));

  p->type = typeCon;

  for ( i = 0; i < NumFunc; i++)
    if ( fun_sym_tbl[i].type == 0 )
      if ( strcmp(fun, fun_sym_tbl[i].name) == 0 )
	{ 
	  funcID = i;
	  break;
	}

  if ( funcID == -1 )
    cdoAbort("function %s not available!", fun);

  p->u.con.value = fun_sym_tbl[funcID].func(p1->u.con.value);

  return (p);
}

static
nodeType *ex_fun_var(char *fun, nodeType *p1)
{
  nodeType *p;
  long ngp, i;
  long nlev;
  int gridID, zaxisID;
  int funcID = -1;
  int nmiss;
  double missval;

  gridID  = p1->gridID;
  zaxisID = p1->zaxisID;
  nmiss   = p1->nmiss;
  missval = p1->missval;

  ngp  = gridInqSize(gridID);
  nlev = zaxisInqSize(zaxisID);

  p = (nodeType *) malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->tmpvar   = 1;
  p->u.var.nm = strdupx("tmp");
  p->gridID   = gridID;
  p->zaxisID  = zaxisID;
  p->missval  = missval;

  p->data = (double *) malloc(ngp*nlev*sizeof(double));

  for ( i = 0; i < NumFunc; i++)
    if ( strcmp(fun, fun_sym_tbl[i].name) == 0 )
      { 
	funcID = i;
	break;
      }

  if ( funcID == -1 )
    cdoAbort("function %s not available!", fun);

  if ( nmiss > 0 )
    {
      for ( i = 0; i < ngp*nlev; i++ )
	{
	  errno = -1;
	  p->data[i] = DBL_IS_EQUAL(p1->data[i], missval) ? missval : fun_sym_tbl[funcID].func(p1->data[i]);
	  if ( errno == EDOM || errno == ERANGE ) p->data[i] = missval;
	  else if ( isnan(p->data[i]) )  p->data[i] = missval;
	}
    }
  else
    {
      for ( i = 0; i < ngp*nlev; i++ )
	{
	  errno = -1;
	  p->data[i] = fun_sym_tbl[funcID].func(p1->data[i]);
	  if ( errno == EDOM || errno == ERANGE ) p->data[i] = missval;
	  else if ( isnan(p->data[i]) )  p->data[i] = missval;
	}
    }

  nmiss = 0;
  for ( i = 0; i < ngp*nlev; i++ )
    if ( DBL_IS_EQUAL(p->data[i], missval) ) nmiss++;

  p->nmiss = nmiss;

  if ( p1->tmpvar ) free(p1->data);

  return (p);
}

static
nodeType *ex_fun(char *fun, nodeType *p1)
{
  nodeType *p = NULL;

  if ( p1->type == typeVar )
    {
      p = ex_fun_var(fun, p1);
      if ( cdoVerbose ) printf("\t%s (%s)\n", fun, p1->u.var.nm);
    }
  else if ( p1->type == typeCon )
    {
      p = ex_fun_con(fun, p1);
      if ( cdoVerbose ) printf("\t%s (%g)\n", fun, p1->u.con.value);
    }
  else
    cdoAbort("Internal problem!");

  return (p);
}

static
nodeType *ex_uminus_var(nodeType *p1)
{
  nodeType *p;
  long ngp, i;
  long nlev;
  int nmiss;
  int gridID, zaxisID;
  double missval;

  gridID   = p1->gridID;
  zaxisID  = p1->zaxisID;
  nmiss    = p1->nmiss;
  missval  = p1->missval;

  ngp  = gridInqSize(gridID);
  nlev = zaxisInqSize(zaxisID);

  p = (nodeType *) malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->tmpvar   = 1;
  p->u.var.nm = strdupx("tmp");
  p->gridID   = gridID;
  p->zaxisID  = zaxisID;
  p->missval  = missval;

  p->data = (double *) malloc(ngp*nlev*sizeof(double));

  if ( nmiss > 0 )
    {
      for ( i = 0; i < ngp*nlev; i++ )
	p->data[i] = DBL_IS_EQUAL(p1->data[i], missval) ? missval : -(p1->data[i]);
    }
  else
    {
      for ( i = 0; i < ngp*nlev; i++ )
	p->data[i] = -(p1->data[i]);
    }

  p->nmiss = nmiss;
  
  return (p);
}

static
nodeType *ex_uminus_con(nodeType *p1)
{
  nodeType *p;

  p = (nodeType *) malloc(sizeof(nodeType));

  p->type = typeCon;

  p->u.con.value = -(p1->u.con.value);

  return (p);
}

static
nodeType *ex_uminus(nodeType *p1)
{
  nodeType *p = NULL;

  if ( p1->type == typeVar )
    {
      p = ex_uminus_var(p1);
      if ( cdoVerbose ) printf("\t- (%s)\n", p1->u.var.nm);
    }
  else if ( p1->type == typeCon )
    {
      p = ex_uminus_con(p1);
      if ( cdoVerbose ) printf("\t- (%g)\n", p1->u.con.value);
    }
  else
    cdoAbort("Internal problem!");

  return (p);
}


int exNode(nodeType *p, parse_parm_t *parse_arg)
{
  int k;              /* child number */

  if ( ! p ) return(0);

  /* node is leaf */
  if ( p->type == typeCon || p->type == typeVar || p->u.opr.nops == 0 )
    {
      return (0);
    }

  /* node has children */
  for ( k = 0; k < p->u.opr.nops; k++ )
    {
      exNode(p->u.opr.op[k], parse_arg);
    }

  return (0);
}


nodeType *expr_run(nodeType *p, parse_parm_t *parse_arg)
{
  int gridID1 = -1, zaxisID1 = -1, tsteptype1 = -1;
  double missval = 0;
  char varname[256];
  int varID, nvars;
  nodeType *rnode = NULL;

  if ( ! p ) return (rnode);

  /*  if ( ! parse_arg->init ) { exNode(p, parse_arg); return (0); } */

  switch ( p->type )
    {
    case typeCon:       
      if ( parse_arg->init )
	{
	  if ( parse_arg->debug )
	    printf("\tpush\t%g\n", p->u.con.value);
	}
      else
	{
	  rnode = p;
	}

      break;
    case typeVar:
      /*    if ( parse_arg->init ) */
	{
	  if ( parse_arg->debug )
	    printf("\tpush\t%s\n", p->u.var.nm);

	  nvars = vlistNvars(parse_arg->vlistID1);
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      vlistInqVarName(parse_arg->vlistID1, varID, varname);
	      if ( strcmp(varname, p->u.var.nm) == 0 ) break;
	    }

	  if ( varID == nvars )
	    {
	      cdoAbort("Variable >%s< not found!", p->u.var.nm);
	    }
	  else
	    {
	      if ( varID >= MAX_VARS ) cdoAbort("Too many parameter (limit=%d)!", MAX_VARS);

	      if ( parse_arg->var_needed[varID] == 0 )
		{

		  parse_arg->var[varID] = strdupx(p->u.var.nm);
		  parse_arg->varID[varID] = varID;
		  parse_arg->var_needed[varID] = 1;
		}

	      gridID1  = vlistInqVarGrid(parse_arg->vlistID1, varID);
	      zaxisID1 = vlistInqVarZaxis(parse_arg->vlistID1, varID);
	      tsteptype1  = vlistInqVarTsteptype(parse_arg->vlistID1, varID);
	      missval  = vlistInqVarMissval(parse_arg->vlistID1, varID);

	      parse_arg->missval2 = missval;

	      if ( parse_arg->gridID2 == -1 )
		parse_arg->gridID2 = gridID1;

	      if ( parse_arg->zaxisID2 == -1 )
		parse_arg->zaxisID2 = zaxisID1;

	      if ( parse_arg->tsteptype2 == -1 || parse_arg->tsteptype2 == TSTEP_CONSTANT )
		parse_arg->tsteptype2 = tsteptype1;
	    }
	}
	/* else */
	{ 
	  if ( parse_arg->debug )
	    printf("%s %d %d %d\n", p->u.var.nm, varID, gridID1, zaxisID1);
	  p->gridID  = gridID1;
	  p->zaxisID = zaxisID1;
	  p->missval = missval;
          p->nmiss   = 0;
	  if ( ! parse_arg->init )
	    {
	      p->data  = parse_arg->vardata1[varID];
	      p->nmiss = parse_arg->nmiss[varID];
	    }
	  p->tmpvar  = 0;
	  rnode = p;
	}

      break;
    case typeFun:
      if ( parse_arg->init )
	{
	  expr_run(p->u.fun.op, parse_arg);

	  if ( parse_arg->debug )
	    printf("\tcall \t%s\n", p->u.fun.name);
	}
      else
	{
	  rnode = ex_fun(p->u.fun.name, expr_run(p->u.fun.op, parse_arg));
	}
      break;
    case typeOpr:
      switch( p->u.opr.oper )
	{
        case '=':
	  parse_arg->gridID2  = -1;
	  parse_arg->zaxisID2 = -1;
          parse_arg->tsteptype2  = -1;

	  rnode = expr_run(p->u.opr.op[1], parse_arg);

	  if ( parse_arg->init )
	    {
	      if ( parse_arg->debug )
		printf("\tpop\t%s\n", p->u.opr.op[0]->u.var.nm);
	      /*
	      if ( p->u.opr.op[1]->type != typeVar )
		cdoAbort("Operand not variable!");
	      */
	      if ( parse_arg->gridID2 == -1 || parse_arg->zaxisID2 == -1 || parse_arg->tsteptype2 == -1 )
		cdoAbort("Operand not variable!");

	      varID = vlistDefVar(parse_arg->vlistID2, parse_arg->gridID2, parse_arg->zaxisID2, parse_arg->tsteptype2);
	      vlistDefVarName(parse_arg->vlistID2, varID, p->u.opr.op[0]->u.var.nm);
	      vlistDefVarMissval(parse_arg->vlistID2, varID, parse_arg->missval2);
	    }
	  else
	    {
	      if ( parse_arg->debug )
		printf("\tpop\t%s\t%s\n", p->u.opr.op[0]->u.var.nm, rnode->u.var.nm);

	      nvars = vlistNvars(parse_arg->vlistID2);
	      for ( varID = 0; varID < nvars; varID++ )
		{
		  vlistInqVarName(parse_arg->vlistID2, varID, varname);
		  if ( strcmp(varname, p->u.opr.op[0]->u.var.nm) == 0 ) break;
		}

	      if ( varID == nvars )
		{
		  cdoAbort("variable >%s< not found!", p->u.opr.op[0]->u.var.nm);
		}
	      else
		{
		  parse_arg->gridID2  = vlistInqVarGrid(parse_arg->vlistID2, varID);
		  parse_arg->zaxisID2 = vlistInqVarZaxis(parse_arg->vlistID2, varID);
		  parse_arg->tsteptype2  = vlistInqVarTsteptype(parse_arg->vlistID2, varID);
		  missval  = vlistInqVarMissval(parse_arg->vlistID2, varID);
	      
		  p->gridID  = parse_arg->gridID2;
		  p->zaxisID = parse_arg->zaxisID2;
		  p->missval = missval;
		  p->data    = parse_arg->vardata2[varID];
		  p->tmpvar  = 0;

		  ex_copy(p, rnode);

		  if ( rnode->tmpvar ) free(rnode->data);
		}
	    }

	  break;
        case UMINUS:    
	  if ( parse_arg->init )
	    {
	      expr_run(p->u.opr.op[0], parse_arg);

	      if ( parse_arg->debug )
		printf("\tneg\n");
	    }
	  else
	    {
	      rnode = ex_uminus(expr_run(p->u.opr.op[0], parse_arg));
	    }

	  break;
        default:
	  if ( parse_arg->init )
	    {
	      expr_run(p->u.opr.op[0], parse_arg);
	      expr_run(p->u.opr.op[1], parse_arg);
	      if ( parse_arg->debug )
		switch( p->u.opr.oper )
		  {
		  case '+':  printf("\tadd\n"); break;
		  case '-':  printf("\tsub\n"); break;
		  case '*':  printf("\tmul\n"); break;
		  case '/':  printf("\tdiv\n"); break;
		  case '<':  printf("\tcompLT\n"); break;
		  case '>':  printf("\tcompGT\n"); break;
		  case GE:   printf("\tcompGE\n"); break;
		  case LE:   printf("\tcompLE\n"); break;
		  case NE:   printf("\tcompNE\n"); break;
		  case EQ:   printf("\tcompEQ\n"); break;
		  }
	    }
	  else
	    {
	      rnode = expr(p->u.opr.oper, expr_run(p->u.opr.op[0], parse_arg),
			                  expr_run(p->u.opr.op[1], parse_arg));
	    }
        }
    }

  return (rnode);
}
