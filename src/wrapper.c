#include <Rinternals.h>
#include <R.h>

int kcluster(double *x, int n, int p, int k, int iter_max,
    double *kernel_matrix, double *mu, double *sse, int *ic1, int heuristic);
void get_kernel_matrix(double *x, int n, int p, double h1, double h2,
                       double  (*kernel)(int, int, double*, int, int, double, double),
                       double *kernel_matrix);
void center_kernel_matrix(double *kernel_matrix, int n);
double kernel_gaussian(int i, int j, double *x, int n, int p, double sigmasq, double blank);
double kernel_poly(int i, int j, double *x, int n, int p, double h, double a);
double kernel_sigmoid(int i, int j, double *x, int n, int p, double theta0, double theta1);

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
SEXP kkmeans(SEXP data, SEXP centers, SEXP kern, SEXP param, SEXP iter_max, SEXP init_centers,
             SEXP method, SEXP kmat)
{
  SEXP dims = getAttrib(data, R_DimSymbol);
  SEXP cluster_out;
  SEXP mu_out;
  SEXP sse_out;
  SEXP ret_list;

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
  int     j;
  int     n         = INTEGER(dims)[0];
  int     p         = INTEGER(dims)[1];
  int    *imax_ptr  = INTEGER(iter_max);
  int    *k_ptr     = INTEGER(centers);
  /* 1 == OT-QT, 2 == MacQueen */
  int    *heur_ptr  = INTEGER(method);
  double *h_ptr     = REAL(param);

  int    imax      = *imax_ptr;
  int    k         = *k_ptr;
  int    heuristic = *heur_ptr;
  double h         = *h_ptr;

  double (*kernel)(int, int, double[], int, int, double, double);
  kernel = NULL;

  const char *kern_string = CHAR(STRING_ELT(kern, 0));

  if (strcmp(kern_string, "gaussian") == 0 || strcmp(kern_string, "g") == 0)
    kernel = &kernel_gaussian;
  else if (strcmp(kern_string, "poly") == 0 || strcmp(kern_string, "p") == 0)
    kernel = &kernel_poly;
  else if (strcmp(kern_string, "sigmoid") == 0 || strcmp(kern_string, "s") == 0)
    kernel = &kernel_sigmoid;
  else
    error("Error: `kern` is not recognized.");

  double *x             = (double *) S_alloc(n * p, sizeof(double));
  double *mu            = (double *) S_alloc(k * p, sizeof(double));
  double *sse           = (double *) S_alloc(k, sizeof(double));
  double *kernel_matrix = (double *) S_alloc(n * n, sizeof(double));
  int    *ic1 = (int *) S_alloc(n, sizeof(int));

  /* set the variables in C */
  x = REAL(data);

  ic1 = INTEGER(init_centers);

  /*  */

  int n_iter;
  SEXP niter;

  /* get_kernel_matrix(x, n, p, h, kernel, kernel_matrix); */
  kernel_matrix = REAL(kmat);
  /* center_kernel_matrix(kernel_matrix, n); */
  /* kcluster(x, n, p, k, imax, kernel_matrix, mu, sse, ic1, heuristic); */
  n_iter = kcluster(x, n, p, k, imax, kernel_matrix, mu, sse, ic1, heuristic);

  /* Create some return values for R */
  PROTECT(cluster_out = allocVector(INTSXP , n));
  PROTECT(mu_out      = allocMatrix(REALSXP, k, p));
  PROTECT(sse_out     = allocVector(REALSXP, k));
  PROTECT(niter       = allocVector(INTSXP , 1));
  PROTECT(ret_list    = allocVector(VECSXP , 4));

  double *pmu;
  double *psse;
  int    *pcluster;

  INTEGER(niter)[0] = n_iter;

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
  SET_VECTOR_ELT(ret_list, 3, niter);


  // TODO: it's easier to set names in R, but it might be better to set them
  // here. Look into it.
  UNPROTECT(5);

  return ret_list;
}
