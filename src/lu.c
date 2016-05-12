#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define POLYBENCH_TIME 1

#include "polybench.h"

static
void init_array (int n,
   double A[ n + 0][n + 0])
{
  int i, j;

  for (i = 0; i < n; i++)
    {
      for (j = 0; j <= i; j++)
 A[i][j] = (double)(-j % n) / n + 1;
      for (j = i+1; j < n; j++) {
 A[i][j] = 0;
      }
      A[i][i] = 1;
    }



  int r,s,t;
  double (*B)[n + 0][n + 0]; B = (double(*)[n + 0][n + 0])polybench_alloc_data ((n + 0) * (n + 0), sizeof(double));;
  for (r = 0; r < n; ++r)
    for (s = 0; s < n; ++s)
      (*B)[r][s] = 0;
  for (t = 0; t < n; ++t)
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
 (*B)[r][s] += A[r][t] * A[s][t];
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
 A[r][s] = (*B)[r][s];
  free((void*)B);;

}




static
void print_array(int n,
   double A[ n + 0][n + 0])

{
  int i, j;

  fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
  fprintf(stderr, "begin dump: %s", "A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf (stderr, "\n");
      fprintf (stderr, "%0.2lf ", A[i][j]);
    }
  fprintf(stderr, "\nend   dump: %s\n", "A");
  fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}




static
void kernel_lu(int n,
        double A[ n + 0][n + 0])
{
  int i, j, k;

#pragma scop
  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
       for (k = 0; k < j; k+=10) {
          A[i][j] -= A[i][k] * A[k][i];
	  A[i][j] -= A[i][k+1] * A[k+1][i];
          A[i][j] -= A[i][k+2] * A[k+2][i];
          A[i][j] -= A[i][k+3] * A[k+3][i];
          A[i][j] -= A[i][k+4] * A[k+4][i];
          A[i][j] -= A[i][k+5] * A[k+5][i];
          A[i][j] -= A[i][k+6] * A[k+6][i];
          A[i][j] -= A[i][k+7] * A[k+7][i];      
 	  A[i][j] -= A[i][k+8] * A[k+8][i];
	  A[i][j] -= A[i][k+9] * A[k+9][i];


 }
        A[i][j] /= A[j][j];
    }
   for (j = i; j < n;j++) {
       for (k = 0; k < i; k+=10) {
          A[i][j] -= A[i][k] * A[k][j];
          A[i][j] -= A[i][k+1] * A[k+1][j];
          A[i][j] -= A[i][k+2] * A[k+2][j];
          A[i][j] -= A[i][k+3] * A[k+3][j];
          A[i][j] -= A[i][k+4] * A[k+4][j];
          A[i][j] -= A[i][k+5] * A[k+5][j];
          A[i][j] -= A[i][k+6] * A[k+6][j];
          A[i][j] -= A[i][k+7] * A[k+7][j];
          A[i][j] -= A[i][k+8] * A[k+8][j];
          A[i][j] -= A[i][k+9] * A[k+9][j];
       }
    }
  }


#pragma endscop
}


int main(int argc, char** argv)
{

  int n = 2000;


  double (*A)[n + 0][n + 0]; A = (double(*)[n + 0][n + 0])polybench_alloc_data ((n + 0) * (n + 0), sizeof(double));;


  init_array (n, *A);


  polybench_timer_start();;


  kernel_lu (n, *A);


  polybench_timer_stop();;
  polybench_timer_print();;



  if (argc > 42 && ! strcmp(argv[0], "")) print_array(n, *A);


  free((void*)A);;

  return 0;
}
