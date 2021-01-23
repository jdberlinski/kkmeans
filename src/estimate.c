#include <Rinternals.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>



void kcluster(double *x, int n, int p, int k, int iter_max,
    double *kernel_matrix, double *mu, double *sse, int *ic1);
void get_kernel_matrix(double *x, int n, int p, double h,
                       double  (*kernel)(int, int, double*, int, int, double),
                       double *kernel_matrix);
double kernel_gaussian(int i, int j, double *x, int n, int p, double sigmasq);
double param_search(double *x, int n, int p, int k, int imax, double *kernel_matrix,
                    double *mu, double *sse, int *ic1, int est_len, double h);

//' An Efficient Kernel K-Means Algorithm
//' @name kkmeans_est
//'
//' @description Performs kernel k-means with the Gaussian kernel using an
//' algorithm similar to Hartigan and Wong's k-means algorithm. The algorithm
//' also estimates the bandwidth parameter sigma.
//'
//' @param data data to cluster
//' @param k the number of clusters
//' @param depth how many iterations for estimation
//' @export
SEXP kkmeans_est(SEXP data, SEXP centers, SEXP depth, SEXP init, SEXP iter_max)
{
  SEXP dims = getAttrib(data, R_DimSymbol);
  SEXP cluster_out;
  SEXP mu_out;
  SEXP sse_out;
  SEXP ret_list;
  SEXP end_sigma;

  if (!isReal(data))
    error("Error: `data` must be a double vector.");
  if (!isInteger(centers))
    error("Error: `k` must be an integer of length one.");
  if (length(centers) != 1)
    error("Error: `k` must be an integer of length one.");
  if (!isInteger(depth))
    error("Error: `depth` must be an integer of length one.");
  if (length(depth) != 1)
    error("Error: `depth` must be an integer of length one.");
  if (!isReal(init))
    error("Error: `init` must be a double of length one.");
  if (length(init) != 1)
    error("Error: `init` must be a double of length one.");
  if (!isInteger(iter_max))
    error("Error: `iter_max` must be an integer of length one.");
  if (length(iter_max) != 1)
    error("Error: `iter_max` must be an integer of length one.");


  int     n         = INTEGER(dims)[0];
  int     p         = INTEGER(dims)[1];
  int    *imax_ptr  = INTEGER(iter_max);
  int    *depth_ptr = INTEGER(depth);
  int    *k_ptr     = INTEGER(centers);
  double *init_ptr  = REAL(init);

  int    imax    = *imax_ptr;
  int    est_len = *depth_ptr;
  int    k       = *k_ptr;
  double h       = *init_ptr;

  double (*kernel)(int, int, double[], int, int, double);
  kernel = &kernel_gaussian;

  double *x             = (double *) S_alloc(n * p, sizeof(double));
  double *mu            = (double *) S_alloc(k * p, sizeof(double));
  double *sse           = (double *) S_alloc(k, sizeof(double));
  double *kernel_matrix = (double *) S_alloc(n * n, sizeof(double));
  int    *ic1 = (int *) S_alloc(n, sizeof(int));

  x = REAL(data);

  /* Do the estimation and clustering */
  get_kernel_matrix(x, n, p, h, kernel, kernel_matrix);

  int converged = 0;
  int *ic1_copy = (int *) S_alloc(n, sizeof(int));

  int j;
  int sum = 0;
  double final_sigma = h;
  while (!converged)
  {

    memcpy(ic1_copy, ic1, n * sizeof(int));
    get_kernel_matrix(x, n, p, final_sigma, kernel, kernel_matrix);
    final_sigma = param_search(x, n, p, k, imax, kernel_matrix, mu, sse, ic1, est_len, final_sigma);

    sum = 0;
    for (j = 0; j < n; j++)
      sum += ic1[j] == ic1_copy[j];

    converged = (sum == n) + (sum == 0);
  }

  kcluster(x, n, p, k, 1, kernel_matrix, mu, sse, ic1);


  PROTECT(cluster_out  = allocVector(INTSXP , n));
  PROTECT(mu_out       = allocMatrix(REALSXP, k, p));
  PROTECT(sse_out      = allocVector(REALSXP, k));
  PROTECT(end_sigma    = allocVector(REALSXP , 1));
  PROTECT(ret_list     = allocVector(VECSXP , 4));

  REAL(end_sigma)[0] = final_sigma;

  double *pmu;
  double *psse;
  int    *pcluster;

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
  SET_VECTOR_ELT(ret_list, 3, end_sigma);

  // TODO: it's easier to set names in R, but it might be better to set them
  // here. Look into it.
  UNPROTECT(5);

  return ret_list;
}

double param_search(double *x,
                  int     n,
                  int     p,
                  int     k,
                  int     imax,
                  double *kernel_matrix,
                  double *mu,
                  double *sse,
                  int    *ic1,
                  int     est_len,
                  double h)
{
  int i, j;
  int sum = 0;

  double sigma = h;
  int numerator = 1;
  int denominator = 1;

  double *k_prime = (double *) S_alloc(n * n, sizeof(double));
  int    *p_prime = (int *) S_alloc(n, sizeof(int));

  kcluster(x, n, p, k, imax, kernel_matrix, mu, sse, ic1);
  /* for (i = 0; i < n; i++) */
  /*   p_prime[i] = ic1[i]; */
  memcpy(p_prime, ic1, n * sizeof(int));
  memcpy(k_prime, kernel_matrix, n*n * sizeof(double));
  /* for (i = 0; i < n*n; i++) */
  /*   k_prime[i] = kernel_matrix[i]; */

  int done = 0;
  int changes = 0;
  int any_changes = 0;

  while (!done)
  {
    numerator *= 2;
    denominator *= 2;

    for (j = 0; j < n*n; j++)
      k_prime[j] = sqrt(k_prime[j]);

    if (!changes)
    {
      /* Rprintf("No changes, reducing sigma\n"); */
      numerator++;
      for (j = 0; j < n*n; j++)
        if (kernel_matrix[j] > 0) kernel_matrix[j] *= k_prime[j];
    }
    else
    {
      /* Rprintf("Found change, increasing sigma (current value: )\n"); */
      numerator--;
      for (j = 0; j < n*n; j++)
        kernel_matrix[j] /= k_prime[j];
      any_changes = 1;
    }

    kcluster(x, n, p, k, imax, kernel_matrix, mu, sse, ic1);
    sum = 0;
    for (j = 0; j < n; j++)
      sum += ic1[j] == p_prime[j];
    changes = (sum == n) + (sum == 0);
    if (changes)
      sigma = h * denominator / numerator;
    memcpy(p_prime, ic1, n * sizeof(int));

    done = log2(denominator) >= est_len;
  }

  if (!any_changes)
    sigma = h * denominator / numerator;
  else if (!changes)
    for (j = 0; j < n*n; j++)
      if (kernel_matrix[j] > 0) kernel_matrix[j] *= k_prime[j];

  return sigma;

}
