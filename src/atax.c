#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define POLYBENCH_TIME 1

#include "polybench.h"

static
void init_array (int m, int n,
   double A[ m + 0][n + 0],
   double x[ n + 0])
{
  int i, j;
  double fn;
  fn = (double)n;

  for (i = 0; i < n; i++)
      x[i] = 1 + (i / fn);
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      A[i][j] = (double) ((i+j) % n) / (5*m);
}




static
void print_array(int n,
   double y[ n + 0])

{
  int i;

  fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
  fprintf(stderr, "begin dump: %s", "y");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (stderr, "\n");
    fprintf (stderr, "%0.2lf ", y[i]);
  }
  fprintf(stderr, "\nend   dump: %s\n", "y");
  fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}




static
void kernel_atax(int m, int n,
   double A[ m + 0][n + 0],
   double x[ n + 0],
   double y[ n + 0],
   double tmp[ m + 0])
{
  int i, j;

#pragma scop
  for (i = 0; i < n; i++)
    y[i] = 0;
  for (i = 0; i < m; i++)
    {
      tmp[i] = 0.0;
      for (j = 0; j < n; j++)
 tmp[i] = tmp[i] + A[i][j] * x[j];
      for (j = 0; j < n; j++)
 y[j] = y[j] + A[i][j] * tmp[i];
    }
#pragma endscop

}


int main(int argc, char** argv)
{

  int m = 1900;
  int n = 2100;


  double (*A)[m + 0][n + 0]; A = (double(*)[m + 0][n + 0])polybench_alloc_data ((m + 0) * (n + 0), sizeof(double));;
  double (*x)[n + 0]; x = (double(*)[n + 0])polybench_alloc_data (n + 0, sizeof(double));;
  double (*y)[n + 0]; y = (double(*)[n + 0])polybench_alloc_data (n + 0, sizeof(double));;
  double (*tmp)[m + 0]; tmp = (double(*)[m + 0])polybench_alloc_data (m + 0, sizeof(double));;


  init_array (m, n, *A, *x);


  polybench_timer_start();


  kernel_atax (m, n,
        *A,
        *x,
        *y,
        *tmp);


  polybench_timer_stop();;
  polybench_timer_print();;



  if (argc > 42 && ! strcmp(argv[0], "")) print_array(n, *y);


  free((void*)A);;
  free((void*)x);;
  free((void*)y);;
  free((void*)tmp);;

  return 0;
}
