#ifndef _MERGE_SORT2_H_
#define _MERGE_SORT2_H_

/* MERGE SORT DEFINES */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cdo.h"

#if defined (_OPENMP)
#include <omp.h>
#endif

static int MERGE_SORT_LIMIT_SIZE = 16384; 
static int first_sort_iter_call = 1;
static double merge_time;

void sort_iter_single(long num_links, double *restrict add1, int parent);

#endif
