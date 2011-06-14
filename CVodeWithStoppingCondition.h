#include <cvode.h>

int CVodeWithStoppingCondition
  (void *cvode_mem, realtype tout, N_Vector yout, 
   realtype *tret, int itask,
   int (*stopping_condition)(N_Vector, realtype, void *));
