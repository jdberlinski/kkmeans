/* @file kernels.c
 * @auth Josh Berlinski
 * @date 2020/05/18
 *
 * File containing all of the kernel functions for the algorithm
 *
 */

#include <math.h>
#include <Rmath.h>

int is_na(double value);
double kernel_gaussian(int i, int j, double *x, int n, int p, double sigmasq, double blank);
double kernel_poly(int i, int j, double *x, int n, int p, double h, double a);
double kernel_sigmoid(int i, int j, double *x, int n, int p, double theta0, double theta1);

int is_na(double value) {
  return value != value;
}

/*
 *  Gaussian RBF kernel
 *
 *  exp( -0.5 * ||x - y||^2 / sigmasq )
 *
 * @param i the index of the first observation
 * @param j index of the second
 * @param x the data vector
 * @param p the dimension
 * @param sigmasq tuning parameter
 * @return the kernel evaluated at the two points
 */
double kernel_gaussian(int i, int j, double *x, int n, int p, double sigmasq, double blank)
{
  double norm = 0.;
  for (int r = 0 ; r < p ; r++)
  {
    if (!is_na(x[i + r*n]) && !is_na(x[j + r*n]))
      norm +=  (x[i + r*n] - x[j + r*n]) * (x[i + r*n] - x[j + r*n]);
  }

  return (exp(-0.5 * norm / sigmasq));
}

/*
 *  Polynomial kernel
 *
 *  ( x'y )^p
 *
 * @param i the index of the first observation
 * @param j index of the second
 * @param x the data vector
 * @param p the dimension
 * @param p tuning parameter
 * @return the kernel evaluated at the two points
 */
double kernel_poly(int i, int j, double *x, int n, int p, double h, double a)
{
  double prod = 0.;
  for (int r = 0; r < p; r++)
  {
    if (!is_na(x[i + r*n]) && !is_na(x[j + r*n]))
      prod += x[i + r*n] * x[j + r*n];
  }

  return (pow(prod + a, h));
}

/*
 * Sigmoid kernel
 *
 * tanh(theta0 * x'y + theta1)
 */
double kernel_sigmoid(int i, int j, double *x, int n, int p, double theta0, double theta1) {
  double prod = 0.;

  for (int r = 0; r < p; r++) {
    if (!is_na(x[i + r*n]) && !is_na(x[j + r*n]))
      prod += x[i + r*n] * x[j + r*n];
  }

  return tanh(theta0 * prod + theta1);
}
