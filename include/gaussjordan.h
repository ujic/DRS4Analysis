#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define FREE_ARG char* 
#define NR_END 1 

void nrerror(char error_text[]);
int *ivector(long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
void free_ivector(int *v, long nl, long nh);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void gaussj(double **a, int n, double **b, int m); // this is added from nr.h
#else /* ANSI */
/* traditional - K&R */
void nrerror();


#endif /* _NR_UTILS_H_ */


