#include <Rinternals.h>
#include <R.h>
void get_kernel_matrix(double *x,
                       int     n,
                       int     p,
                       double  h,
                       double  (*kernel)(int, int, double*, int, int, double),
                       double *kernel_matrix);
double kernel_gaussian(int i, int j, double *x, int n, int p, double sigmasq);
double kernel_poly(int i, int j, double *x, int n, int p, double h);

//' Get the kernel matrix for a dataset
//' @name get_k_matrix
//'
//' @description Given a dataset, kernel function, and tuning parameter, will
//' return the n x n kernel matrix
//'
//' @param data data vector
//' @param kern the kernel to use, one of ('gaussian', 'poly'), can use first
//' letter
//' @param params parameters to pass to kernel function.
//' @export
SEXP get_k_matrix(SEXP data, SEXP kern, SEXP param)
{
  SEXP dims = getAttrib(data, R_DimSymbol);
  SEXP kmatrix;

  if (!isReal(data))
    error("Error: `data` must be a double vector.");
  if (!isString(kern) || length(kern) != 1)
    error("Error: `kern` must be a single string");
  if (!isReal(param))
    error("Error: `param` must be a double of length one.");
  if (length(param) != 1)
    error("Error: `param` must be a double of length one.");

  /* get the variables from R */
  int     j;
  int     n        = INTEGER(dims)[0];
  int     p        = INTEGER(dims)[1];
  double *h_ptr    = REAL(param);

  double h    = *h_ptr;

  double (*kernel)(int, int, double[], int, int, double);
  kernel = NULL;

  const char *kern_string = CHAR(STRING_ELT(kern, 0));

  if (strcmp(kern_string, "gaussian") == 0 || strcmp(kern_string, "g") == 0)
    kernel = &kernel_gaussian;
  else if (strcmp(kern_string, "poly") == 0 || strcmp(kern_string, "p") == 0)
    kernel = &kernel_poly;
  else
    error("Error: `kern` is not recognized.");

  double *x   = (double *) S_alloc(n * p, sizeof(double));

  x = REAL(data);

  double *kernel_matrix = (double *) S_alloc(n * n, sizeof(double));

  get_kernel_matrix(x, n, p, h, kernel, kernel_matrix);

  PROTECT(kmatrix = allocMatrix(REALSXP, n, n));

  double *pkmatrix;

  pkmatrix = REAL(kmatrix);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      pkmatrix[i + j*n] = kernel_matrix[i + j*n];

  UNPROTECT(1);

  return kmatrix;
}
