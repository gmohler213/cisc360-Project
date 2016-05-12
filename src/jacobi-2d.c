#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define POLYBENCH_TIME 1
#define MIN4(x, y) (((x+4) < (y)) ? (x+4) : (y))
#define MIN8(x, y) (((x+8) < (y)) ? (x+8) : (y))

#include "polybench.h"

static
void init_array (int n,
   double A[ n + 0][n + 0],
   double B[ n + 0][n + 0])
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
 A[i][j] = ((double) i*(j+2) + 2) / n;
 B[i][j] = ((double) i*(j+3) + 3) / n;
      }
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
      if ((i * n + j) % 20 == 0) fprintf(stderr, "\n");
      fprintf(stderr, "%0.2lf ", A[i][j]);
    }
  fprintf(stderr, "\nend   dump: %s\n", "A");
  fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}




static
void kernel_jacobi_2d(int tsteps,
       int n,
       double A[n][n],
       double B[n][n])

{
  int t, i, j;

  #pragma scop
  for (t = 0; t < tsteps; t++)
    {
      for (i = 1; i < n - 1; i++)
	for (j = 1; j < n - 1; j+=5){
	  B[i][j] = 0.2 * (A[i][j] + A[i][j-1] + A[i][1+j] + A[1+i][j] + A[i-1][j]);
	  B[i][j+1] = 0.2 * (A[i][j+1] + A[i][j] + A[i][2+j] + A[1+i][j+1] + A[i-1][j+1]);
	  B[i][j+2] = 0.2 * (A[i][j+2] + A[i][j+1] + A[i][j+3] + A[1+i][j+2] + A[i-1][j+2]);
	  B[i][j+3] = 0.2 * (A[i][j+3] + A[i][j+2] + A[i][j+4] + A[1+i][j+3] + A[i-1][j+3]);
	  B[i][j+4] = 0.2 * (A[i][j+4] + A[i][j+3] + A[i][j+5] + A[1+i][j+4] + A[i-1][j+4]);
	}
      for (i = 1; i < n - 1; i++)
	for (j = 1; j < n - 1; j+=5){
	  A[i][j] = 0.2 * (B[i][j] + B[i][j-1] + B[i][1+j] + B[1+i][j] + B[i-1][j]);
	  A[i][j+1] = 0.2 * (B[i][j+1] + B[i][j] + B[i][2+j] + B[1+i][j+1] + B[i-1][j+1]);
	  A[i][j+2] = 0.2 * (B[i][j+2] + B[i][j+1] + B[i][j+3] + B[1+i][j+2] + B[i-1][j+2]);
	  A[i][j+3] = 0.2 * (B[i][j+3] + B[i][j+2] + B[i][j+4] + B[1+i][j+3] + B[i-1][j+3]);
	  A[i][j+4] = 0.2 * (B[i][j+4] + B[i][j+3] + B[i][j+5] + B[1+i][j+4] + B[i-1][j+4]);
	}
    }
  #pragma endscop

} 
/*
{
  int t, i, j;

#pragma scop
  for (t = 0; t < tsteps; t++)
    {
      for (i = 1; i < n - 1; i++)
 for (j = 1; j < n - 1; j++)
   B[i][j] = 0.2 * (A[i][j] + A[i][j-1] + A[i][1+j] + A[1+i][j] + A[i-1][j]);
      for (i = 1; i < n - 1; i++)
 for (j = 1; j < n - 1; j++)
   A[i][j] = 0.2 * (B[i][j] + B[i][j-1] + B[i][1+j] + B[1+i][j] + B[i-1][j]);
    }
#pragma endscop

}
*/


int main(int argc, char** argv)
{

  int n = 1300;
  int tsteps = 500;


  double (*A)[n + 0][n + 0]; A = (double(*)[n + 0][n + 0])polybench_alloc_data ((n + 0) * (n + 0), sizeof(double));;
  double (*B)[n + 0][n + 0]; B = (double(*)[n + 0][n + 0])polybench_alloc_data ((n + 0) * (n + 0), sizeof(double));;



  init_array (n, *A, *B);


  polybench_timer_start();;


  kernel_jacobi_2d(tsteps, n, *A, *B);


  polybench_timer_stop();;
  polybench_timer_print();;



  if (argc > 42 && ! strcmp(argv[0], "")) print_array(n, *A);


  free((void*)A);;
  free((void*)B);;

  return 0;
}
