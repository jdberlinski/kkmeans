#include <Rcpp.h>
#include <cstring>

using namespace Rcpp;

void kcluster(double x[], int n, int p, int k, double h, int iter_max,
    double (*kernel)(int, int, double[], int, int, double),
    double mu[], double sse[], int ic1[]);
double kernel_gaussian(int i, int j, double x[], int n, int p, double sigmasq);
double kernel_poly(int i, int j, double x[], int n, int p, double h);

//' @export
// [[Rcpp::export]]
List kkmeans(NumericMatrix data, int k, String kern, double param, int iter_max = 1000) {

  int n = data.nrow();
  int p = data.ncol();

  double* x = (double *) S_alloc(n * p, sizeof(double));
  double* mu = (double *) S_alloc(k * p, sizeof(double));
  double* sse = (double *) S_alloc(k, sizeof(double));
  int* ic1 = (int *) S_alloc(n, sizeof(int));


  x = data.begin();

  double (*kernel)(int, int, double[], int, int, double);

  const char* kern_string = kern.get_cstring();

  if (strcmp(kern_string, "gaussian") == 0 || strcmp(kern_string, "g") == 0)
    kernel = &kernel_gaussian;
  else if (strcmp(kern_string, "poly") == 0 || strcmp(kern_string, "p") == 0)
    kernel = &kernel_poly;
  else
    stop("Error: `kern` is not recognized.");

  kcluster(x, n, p, k, param, iter_max, kernel, mu, sse, ic1);

  // This is a hack because there is no clear way to convert a pointer to a
  // Vector object?
  IntegerVector cluster_out(n);
  for (int i = 0; i < n; i++)
    cluster_out[i] = ic1[i] + 1;

  NumericMatrix mu_out(k, p, mu);

  NumericVector sse_out(k);
  for (int i = 0; i < k; i++)
    sse_out[i] = sse[i];

  List out = List::create(_["cluster"] = cluster_out, _["center"] = mu_out,
      _["sse"] = sse_out);
  return out;
}
