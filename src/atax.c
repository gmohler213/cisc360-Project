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
   double A[ m ][n ],
   double x[ n ],
   double y[ n ],
   double tmp[ m ])
{
  int i, j;

#pragma scop
  for (i = 0; i < n; i+=10){
    y[i] = 0;
    y[i+1] = 0;
    y[i+2] = 0;
    y[i+3] = 0;
    y[i+4] = 0;
    y[i+5] = 0;
    y[i+6] = 0;
    y[i+7] = 0;
    y[i+8] = 0;
    y[i+9] = 0;
  }
  for (i = 0; i < m; i++)
    {
      tmp[i] = 0.0;
      for (j = 0; j < n; j+=10){
	tmp[i] += A[i][j] * x[j];
	tmp[i] += A[i][j+1] * x[j+1];
	tmp[i] += A[i][j+2] * x[j+2];
	tmp[i] += A[i][j+3] * x[j+3];
	tmp[i] += A[i][j+4] * x[j+4];
	tmp[i] += A[i][j+5] * x[j+5];
	tmp[i] += A[i][j+6] * x[j+6];
	tmp[i] += A[i][j+7] * x[j+7];
	tmp[i] += A[i][j+8] * x[j+8];
	tmp[i] += A[i][j+9] * x[j+9];
      }
      for (j = 0; j < n; j+=10){
	y[j] += A[i][j] * tmp[i];
	y[j+1] += A[i][j+1] * tmp[i];
	y[j+2] += A[i][j+2] * tmp[i];
	y[j+3] += A[i][j+3] * tmp[i];
	y[j+4] += A[i][j+4] * tmp[i];
	y[j+5] += A[i][j+5] * tmp[i];
	y[j+6] += A[i][j+6] * tmp[i];
	y[j+7] += A[i][j+7] * tmp[i];
	y[j+8] += A[i][j+8] * tmp[i];
	y[j+9] += A[i][j+9] * tmp[i];
      }
    }

  
#pragma endscop
}
/*
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
*/

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
