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

//' @export
SEXP classify_kkmeans(SEXP test_data, SEXP train_data, SEXP labels, SEXP kern, SEXP param, SEXP param2,
    SEXP ncluster) {
  SEXP dims_train = getAttrib(train_data, R_DimSymbol);
  SEXP dims_test = getAttrib(test_data, R_DimSymbol);

  if (!isReal(test_data))
    error("Error: `test_data` must be a double vector.");
  if (!isReal(train_data))
    error("Error: `train_data` must be a double vector.");
  if (!isString(kern) || length(kern) != 1)
    error("Error: `kern` must be a single string");
  if (!isReal(param))
    error("Error: `param` must be a double of length one.");
  if (length(param) != 1)
    error("Error: `param` must be a double of length one.");
  if (!isReal(param2))
    error("Error: `param` must be a double of length one.");
  if (length(param2) != 1)
    error("Error: `param` must be a double of length one.");

  int n_train = INTEGER(dims_train)[0];
  int p_train = INTEGER(dims_train)[1];
  int n_test = INTEGER(dims_test)[0];
  int p_test = INTEGER(dims_test)[1];
  int p = 0;

  if (p_train != p_test)
    error("Error: Test data must have the number of columns as training data.");
  else
    p = p_train;

  double *h_ptr  = REAL(param);
  double *h2_ptr = REAL(param2);
  double h  = *h_ptr;
  double h2 = *h2_ptr;

  int *k_ptr = INTEGER(ncluster);
  int k = *k_ptr;

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

  double *x = (double *) S_alloc(n_train * p, sizeof(double));
  double *test = (double *) S_alloc(n_test * p, sizeof(double));
  x = REAL(train_data);
  test = REAL(test_data);

  /* we will use an n+1 matrix and add observations to the begining of the vector */
  double *work = (double *) S_alloc((n_train + 1) * p_train, sizeof(double));
  /* the matrix is actually in column-major order so we cannot memcpy */
  for (int i = 1; i <= n_train; i++)
    for (int j = 0; j < p; j++)
      work[i + j*(n_train + 1)] = x[i - 1 + j*n_train];


  int *train_lab = (int *) S_alloc(n_train, sizeof(int));
  train_lab = INTEGER(labels);

 /* kernel matrix for the training data */
  double *kernel_matrix = (double *) S_alloc(n_train * n_train, sizeof(double));
  get_kernel_matrix(x, n_train, p_train, h, h2, kernel, kernel_matrix);

  /* compute cross-kernel for each cluster */
  double *kern_cross = (double *) S_alloc(k, sizeof(double));
  int *n_k = (int *) S_alloc(k, sizeof(int));

  for (int i = 0; i < n_train; i++) {
    n_k[train_lab[i]]++;
    for (int j = 0; j < n_train; j++)
      if (train_lab[i] == train_lab[j])
        kern_cross[train_lab[i]] += kernel_matrix[i + j*n_train];
  }

  for (int i = 0; i < k; i++)
    kern_cross[i] /= (double) (n_k[i] * n_k[i]);

  double self_kern;

  double big = 1.0E+10;
  double current_min;
  double dist;

  int current_cluster = -1;

  double *first_order = (double *) S_alloc(k, sizeof(double));

  int *test_labels = (int *) S_alloc(n_test, sizeof(int));
  /* for each new observation, find the cluster with the lowest distance */
  for (int i = 0; i < n_test; i++) {
    memset(&first_order[0], 0, k * sizeof(double));

    for (int j = 0; j < p; j++)
      work[j*(n_train + 1)] = test[i + j*n_test];

    self_kern = kernel(0, 0, work, n_train + 1, p, h, h2);

    for (int j = 0; j < n_train; j++)
      first_order[train_lab[j]] += kernel(0, j + 1, work, n_train + 1, p, h, h2);

    current_min = big;
    for (int j = 0; j < k; j++) {
      dist = self_kern - 2.0 / n_k[j] * first_order[j] + kern_cross[j];
      if (dist < current_min) {
        current_min = dist;
        current_cluster = j;
      }
    }

    test_labels[i] = current_cluster;
  }

  SEXP cluster_out;
  PROTECT(cluster_out = allocVector(INTSXP, n_test));
  int *pcluster;
  pcluster = INTEGER(cluster_out);
  for  (int i = 0; i < n_test; i++)
    pcluster[i] = test_labels[i] + 1;
  UNPROTECT(1);

  return cluster_out;
}
