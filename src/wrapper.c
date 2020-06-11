#include <Rinternals.h>
#include <R.h>

void kcluster(double *x, int n, int p, int k, double h, int iter_max,
    double (*kernel)(int, int, double*, int, int, double),
    double *mu, double *sse, int *ic1);
double kernel_gaussian(int i, int j, double *x, int n, int p, double sigmasq);
double kernel_poly(int i, int j, double *x, int n, int p, double h);

//' An Efficient Kernel K-Means Algorithm
//' @name kkmeans
//'
//' @description Performs kernel k-means with the specified kernel using an
//' algorithm similar to Hartigan and Wong's k-means algorithm.
//'
//' @param data data to cluster
//' @param k the number of clusters
//' @param kern the kernel to use, one of ('gaussian', 'poly'), can use first
//' letter
//' @param params parameters to pass to kernel function.
//' @export
SEXP kkmeans(SEXP data, SEXP centers, SEXP kern, SEXP param, SEXP iter_max)
{
  SEXP dims = getAttrib(data, R_DimSymbol);
  /* TODO:  */
  /* SEXP param_names = getAttrib(param, R_NamesSymbol); */
  SEXP cluster_out, mu_out, sse_out, ret_list;

  if (!isReal(data))
    error("Error: `data` must be a double vector.");
  if (!isString(kern) || length(kern) != 1)
    error("Error: `kern` must be a single string");
  if (!isInteger(centers))
    error("Error: `k` must be an integer of length one.");
  if (length(centers) != 1)
    error("Error: `k` must be an integer of length one.");
  if (!isReal(param))
    error("Error: `param` must be a double of length one.");
  if (length(param) != 1)
    error("Error: `param` must be a double of length one.");
  if (!isInteger(iter_max))
    error("Error: `iter_max` must be an integer of length one.");
  if (length(iter_max) != 1)
    error("Error: `iter_max` must be an integer of length one.");

  /* get the variables from R */
  int j;
  int n         = INTEGER(dims)[0];
  int p         = INTEGER(dims)[1];
  int *imax_ptr = INTEGER(iter_max);
  int *k_ptr    = INTEGER(centers);
  double *h_ptr = REAL(param);

  int imax = *imax_ptr;
  int k = *k_ptr;
  double h = *h_ptr;

  double (*kernel)(int, int, double[], int, int, double);
  kernel = NULL;

  const char *kern_string = CHAR(STRING_ELT(kern, 0));
  /* const char* sigmasq; */

  if (strcmp(kern_string, "gaussian") == 0 || strcmp(kern_string, "g") == 0)
    kernel = &kernel_gaussian;
  else if (strcmp(kern_string, "poly") == 0 || strcmp(kern_string, "p") == 0)
    kernel = &kernel_poly;
  else
    error("Error: `kern` is not recognized.");

  double *x   = (double *) S_alloc(n * p, sizeof(double));
  double *mu  = (double *) S_alloc(k * p, sizeof(double));
  double *sse = (double *) S_alloc(k, sizeof(double));
  int *ic1    = (int *) S_alloc(n, sizeof(int));

  /* set the variables in C */
  x = REAL(data);

  kcluster(x, n, p, k, h, imax, kernel, mu, sse, ic1);

  /* Create some return values for R */
  PROTECT(cluster_out = allocVector(INTSXP, n) );
  PROTECT(mu_out = allocMatrix(REALSXP, k, p) );
  PROTECT(sse_out = allocVector(REALSXP, k) );
  PROTECT(ret_list = allocVector(VECSXP, 3) );

  double *pmu, *psse;
  int *pcluster;

  pcluster = INTEGER(cluster_out);
  for (j = 0; j < n; j++)
    pcluster[j] = ic1[j] + 1;

  pmu = REAL(mu_out);
  for (j = 0; j < k; j++)
    for (int m = 0; m < p; m++)
      pmu[j + m*k] = mu[j + m*k];

  psse = REAL(sse_out);
  for (j = 0; j < k; j++)
    psse[j] = sse[j];

  SET_VECTOR_ELT(ret_list, 0, cluster_out);
  SET_VECTOR_ELT(ret_list, 1, mu_out);
  SET_VECTOR_ELT(ret_list, 2, sse_out);

  // TODO: it's easier to set names in R, but it might be better to set them
  // here. Look into it.
  UNPROTECT(4);

  return ret_list;
}
