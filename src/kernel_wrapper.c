#include <Rinternals.h>
#include <R.h>
void get_kernel_matrix(double *x,
                       int     n,
                       int     p,
                       double  h,
                       double  h2,
                       double  (*kernel)(int, int, double*, int, int, double, double),
                       double *kernel_matrix);
double kernel_gaussian(int i, int j, double *x, int n, int p, double sigmasq, double blank);
double kernel_laplacian(int i, int j, double *x, int n, int p, double sigmasq, double blank);
double kernel_poly(int i, int j, double *x, int n, int p, double h, double a);
double kernel_sigmoid(int i, int j, double *x, int n, int p, double theta0, double theta1);

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
SEXP get_k_matrix(SEXP data, SEXP kern, SEXP param, SEXP param2)
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
  double *h2_ptr    = REAL(param2);

  double h    = *h_ptr;
  double h2    = *h2_ptr;

  double (*kernel)(int, int, double[], int, int, double, double);
  kernel = NULL;

  const char *kern_string = CHAR(STRING_ELT(kern, 0));

  if (strcmp(kern_string, "gaussian") == 0 || strcmp(kern_string, "g") == 0)
    kernel = &kernel_gaussian;
  else if (strcmp(kern_string, "laplacian") == 0 || strcmp(kern_string, "l") == 0)
    kernel = &kernel_laplacian;
  else if (strcmp(kern_string, "poly") == 0 || strcmp(kern_string, "p") == 0)
    kernel = &kernel_poly;
  else if (strcmp(kern_string, "sigmoid") == 0 || strcmp(kern_string, "s") == 0)
    kernel = &kernel_sigmoid;
  else
    error("Error: `kern` is not recognized.");

  double *x   = (double *) S_alloc(n * p, sizeof(double));

  x = REAL(data);

  double *kernel_matrix = (double *) S_alloc(n * n, sizeof(double));

  // TODO: center kernel matrix
  get_kernel_matrix(x, n, p, h, h2, kernel, kernel_matrix);

  PROTECT(kmatrix = allocMatrix(REALSXP, n, n));

  double *pkmatrix;

  pkmatrix = REAL(kmatrix);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      pkmatrix[i + j*n] = kernel_matrix[i + j*n];

  UNPROTECT(1);

  return kmatrix;
}
