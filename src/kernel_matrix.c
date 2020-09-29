#include <stdlib.h>
#include "kcluster.h"

void get_kernel_matrix(double *x,
                       int     n,
                       int     p,
                       double  h,
                       double  (*kernel)(int, int, double*, int, int, double),
                       double *kernel_matrix)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j <= i; j++)
    {
      kernel_matrix[get_index(i, j, n)] = kernel(i, j, x, n, p, h);
      kernel_matrix[get_index(j, i, n)] = kernel_matrix[get_index(i, j, n)];
    }
  }
  return;
}
