/*
 * a little utility functions from Numerical Recipes
 */
#include "util.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define NR_END 1
#define FREE_ARG char *
#if 0
/* allocate a vector that can be subscripted [nl .. nh] */
double *dvector(long nl, long nh)
{
  double *v = (double *)malloc((size_t)((nh-nl+1+NR_END)*(sizeof(*v))));
  if (!v) perror("allocation failure in vector()");
  return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
{
  free((FREE_ARG)(v+nl-NR_END));
}
#endif

void *memdup(void *ptr, size_t sz)
{
  void *dp = malloc(sz);
  return memcpy(dp, ptr, sz);
}

const char *glob_string(const char *string, const char *sep)
{
  static char *glob_holding;
  static int  holding_no;
  glob_t globbuf;
  glob(string, GLOB_NOSORT | GLOB_NOCHECK, NULL, &globbuf);
  int ccount = 0;
  int sc = strlen(sep);
  for (unsigned int i = 0; i < globbuf.gl_pathc; i++)
    ccount += strlen(globbuf.gl_pathv[i]) + sc;
  if ( ccount > holding_no )
    glob_holding = (char *)realloc(glob_holding, (holding_no = ccount));
  char *cx = glob_holding;
  for (unsigned int i = 0; i < globbuf.gl_pathc; i++)
  {
    cx = strcpy(cx, globbuf.gl_pathv[i]);
    cx = strchr(cx, '\0');
    if ( i < globbuf.gl_pathc - 1 )
    {
      cx = strcpy(cx, sep);
      cx = strchr(cx, '\0');
    }
  }
  return glob_holding;
}
