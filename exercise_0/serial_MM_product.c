/******************************************************************************
* FILE: mm.c
* DESCRIPTION:  
*   This program calculates the product of matrix a[nra][nca] and b[nca][ncb],
*   the result is stored in matrix c[nra][ncb].
*   The max dimension of the matrix is constraint with static array declaration,
*   for a larger matrix you may consider dynamic allocation of the arrays, 
*   but it makes a parallel code much more complicated (think of communication),
*   so this is only optional.
*   
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main (int argc, char *argv[])
{
int     i, j, k;
int N = 1000;
if (argc > 1) N = atoi(argv[1]);
/* for simplicity, set NRA=NCA=NCB=N  */
int NRA=N;
int NCA=N;
int NCB=N;

double  a[NRA][NCA],           /* matrix A to be multiplied */
        b[NCA][NCB],           /* matrix B to be multiplied */
        c[NRA][NCB];           /* result matrix C */


  /*** Initialize matrices ***/
  
  for (i=0; i<NRA; i++)
    for (j=0; j<NCA; j++)
      a[i][j]= i+j;
  
  for (i=0; i<NCA; i++)
    for (j=0; j<NCB; j++)
      b[i][j]= i*j;
  
  for (i=0; i<NRA; i++)
    for (j=0; j<NCB; j++)
      c[i][j]= 0;

  struct timespec t0, t1;
  clock_gettime(CLOCK_MONOTONIC, &t0);

  for (i=0; i<NRA; i++)
    {
    for(j=0; j<NCB; j++)
      for (k=0; k<NCA; k++)
        c[i][j] += a[i][k] * b[k][j];
    }

  clock_gettime(CLOCK_MONOTONIC, &t1);
  double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
  printf("N=%d, time=%f seconds\n", N, elapsed);

  return 0;
}


