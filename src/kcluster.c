/**
 * @file kkmeans.c
 * @auth Josh Berlinski
 *
 *
 * Utility functions to carry out a Hartigan-Wong type algorithm for kernel K
 * means
 */
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <R.h>
#include <stdio.h>
#include <stdlib.h>

#include "kcluster.h"

/*
 * The main clustering algorithm, with the goal of dividing n observations into
 * k clusters such that the WSS (in feature space) is minimized
 *
 * @param *x        - n*p length array containing col-major data
 * @param  n        - number of data points
 * @param  p        - dimension of data
 * @param  k        - number of clusters
 * @param  iter_max - the maximum number of passes allowed
 * @param  kernel_matrix   - the matrix containing all k(x, y) values
 * @param *mu       - k*p array of cluster centers
 * @param *sse      - k array of SSE for each cluster
 * @param *ic1      - n array of closest cluster indicator
 * @param  heuristic - integer indication of which algorithm to run
 */
int kcluster(double *x,
              int     n,
              int     p,
              int     k,
              int     iter_max,
              double *kernel_matrix,
              double *mu,
              double *sse,
              int    *ic1,
              int     heuristic)
{
  /* n_minus is n_k / (n_k - 1),  n_plus is n_k / (n_k + 1)*/
  double *n_minus = (double *) S_alloc(k, sizeof(double));
  double *n_plus  = (double *) S_alloc(k, sizeof(double));

  /* loss is the current contribution to the SSE of an observation  */
  double *loss = (double *) S_alloc(n, sizeof(double));

  /* foi is the sum of the dot products in feature space for an observation with
   * all other observations in its ith closest cluster*/
  double *fo1 = (double *) S_alloc(n, sizeof(double));
  double *fo2 = (double *) S_alloc(n, sizeof(double));

  /* kern_cross is the sum of all pairs of dot products in each cluster */
  double *kern_cross = (double *) S_alloc(k, sizeof(double));

  /* ic2 is the second closest cluster */
  int *ic2 = (int *) S_alloc(n, sizeof(int));

  /* change tracks the last step at which the cluster was updated
   * if in quick transfer stage, add n (since it comes after the optimal
   * transfer stage)  */
  int *change = (int *) S_alloc(k, sizeof(int));

  /* itran is an indicator of a cluster being changed in the quick transfer
   * stage */
  int *itran = (int *) S_alloc(k, sizeof(int));

  /* live is the indicator if a cluster is in the live set */
  /* in the optimal transfer stage, a cluster is in the live set if live > i,
   * since that means it HAS been updated in the past n optimal transfer stages */
  int *live = (int *) S_alloc(k, sizeof(int));

  /* n_k is the number of observations in a cluster */
  int *n_k = (int *) S_alloc(k, sizeof(int));

  /* kernel_matrix stores the value of the kernel function for each pair of
   * observations */
  double big = 1.0E+10;

  double da;
  int ii;

  double fo;


  /* iteration variables */
  /* Generally, i and j correspond to observations while l corresponds
   * to clusters */
  int i, j, l;

  /* n_transfer is the number of transfers that took place in the optimal
   * transfer stage */
  int n_transfer = -1;

  if (iter_max == 1111)
  {
    init_centers(n, k, ic1, ic2, n_k, kernel_matrix);
  }
  else {
    /* First, randomly assign each point to a cluster */
    for (i = 0; i < n; i++)
    {
      n_k[ic1[i]]++;
      ic2[i] = (ic1[i] + 1) % k;
    }
  }

  for (l = 0; l < k; l++)
  {
    kern_cross[l] = 0.;
    if (!n_k[l])
    {
      i = rand_dunif(n);
      ic1[i] = l;
      n_k[l]++;
      ic2[i] = (ic1[i] + 1) % k;
    }
    for (i = 0; i < n; i++)
    {
      if (ic1[i] == l)
      {
        for (j = 0; j < n; j++)
        {
          if (ic1[j] == l)
            kern_cross[l] += kernel_matrix[get_index(i, j, n)];
        }
      }
    }
  }

  for (l = 0; l < k; l++)
  {
    if (!n_k[l])
      error("Error: Cluster %d has 0 observations. Exiting...\n Possible bad starting seed.\n", l);

    if (n_k[l] == 1)
      n_minus[l] = big;
    else
    {
      n_minus[l] = n_k[l];
      n_minus[l] /= (n_k[l] - 1);
    }

    n_plus[l] = (double) n_k[l];
    n_plus[l] /= (n_k[l] + 1);

    itran[l] = 1;
    change[l] = -1;
  }

  int nqtran, notran;
  nqtran = 0;
  notran = 0;
  int npass = 0;
  for (int iter = 0; iter < iter_max; iter++)
  {

    if (heuristic == 1)
    {
      /* The optimal transfer stage passes through the data once, reallocating
       * each point to the cluster that will yield the greatest reduction of
       * within-cluster sum of squares */

      n_transfer = optimal_transfer(x, mu, sse, n_minus, n_plus, n_k, n, p, k,
          ic1, ic2, change, itran, loss, live, kern_cross, kernel_matrix, fo1, fo2);
      notran++;

      /* if no transfer took place in the optimal transfer stage, then stop */
      /* (because none of the points will be in the live set) */

      npass++;
      if (!n_transfer)
        break;

      quick_transfer(x, mu, sse, n_minus, n_plus, n_k, n, p, k, ic1, ic2, change,
          itran, loss, live, kern_cross, kernel_matrix, fo1, fo2);
      nqtran++;

      /* if there are only two clusters, there is no need to re-enter the optimal
       * transfer stage */
      /* if (k == 2) */
      /*   break; */

      /* change needs to be set back to zero before re-entering the optimal
       * transfer stage */
      memset(&change[0], 0, k * sizeof(int));
      /* for (l = 0; l < k; l++) */
      /*   change[l] = 0; */

    /* Rprintf("Number of optimal transfer iterations: %d \n", notran); */
    /* Rprintf("Number of quick transfer iterations: %d \n", nqtran); */
    }
    else if (heuristic == 2)
    {
      n_transfer = macqueen_step(x, mu, sse, n_minus, n_plus, n_k, n, p, k, ic1,
          loss, kern_cross, kernel_matrix, fo1);

      npass++;
      if (!n_transfer) break;
    }
    else if (heuristic == 3)
    {
      n_transfer = lloyd_step(x, mu, sse, n_minus, n_plus, n_k, n, p, k, ic1,
          loss, kern_cross, kernel_matrix, fo1);

      /* now need to update the cluster-specific values */
      for (l = 0; l < k; l++)
      {
        kern_cross[l] = 0.;
        n_k[l] = 0;

        for (i = 0; i < n; i++)
        {
          if (ic1[i] == l)
          {
            n_k[l]++;
            for (j = 0; j < n; j++)
            {
              if (ic1[j] == l)
                kern_cross[l] += kernel_matrix[get_index(i, j, n)];
            }
          }
        }

        if (n_k[l] <= 1)
          n_minus[l] = big;
        else
        {
          n_minus[l] = (double) n_k[l];
          n_minus[l] /= (n_k[l] - 1);
        }

        n_plus[l] = (double) n_k[l];
        n_plus[l] /= (n_k[l] + 1);
      }


      npass++;
      if (!n_transfer) break;
    }
    else if (heuristic == 4)
    {
      /* The optimal transfer stage passes through the data once, reallocating
       * each point to the cluster that will yield the greatest reduction of
       * within-cluster sum of squares */

      n_transfer = optimal_transfer(x, mu, sse, n_minus, n_plus, n_k, n, p, k,
          ic1, ic2, change, itran, loss, live, kern_cross, kernel_matrix, fo1, fo2);
      notran++;

      /* if no transfer took place in the optimal transfer stage, then stop */
      /* (because none of the points will be in the live set) */

      // TODO: change the exit criteria.
      // could use minimum number of transfers
      // could use epsilon cutoff for wss
      npass++;
      if (!n_transfer)
        break;

      /* change needs to be set back to zero before re-entering the optimal
       * transfer stage */
      memset(&change[0], 0, k * sizeof(int));

    }
    else
    {
      error("Incorrect heuristic");
    }
  }

  // TODO: Something here needs to be fixed. When testing the bullseye dataset with sigma = 40ish,
  // the wss returned is negative.
  // I think kern_cross may need to be updated.
  for (l = 0; l < k; l++)
  {
    kern_cross[l] = 0;
    for (i = 0; i < n; i++)
      if (ic1[i] == l)
        for (j = 0; j < n; j++)
          if (ic1[j] == l)
            kern_cross[l] += kernel_matrix[get_index(i, j, n)];
  }
  for (l = 0; l < k; l++)
  {
    sse[l] = 0;
    /* n_k[l] = 0; */
    for (i = 0; i < n; i++)
    {
      fo1[i] = 0;
      if (ic1[i] != l) continue;
      /* n_k[l]++; */
      sse[l] += kernel_matrix[get_index(i, i, n)];
      sse[l] += 1. / (n_k[l] * n_k[l]) * kern_cross[l];
      fo = 0.;
      for (j = 0; j < n; j++)
      {
        if(ic1[j] != l || j == i) continue;
        fo += kernel_matrix[get_index(i, j, n)];
      }

      fo1[i] = fo;
      fo /= n_k[l];
      sse[l] -= 2*fo;
    }
    /* for (j = 0; j < n; j++) */
    /*   sse[l] += ((double) n_k[l]) * fo1[j]; */
  }

  // TODO: change to memset (faster)
  for (l = 0; l < k; l++)
  {
    /* sse[l] = 0.; */
    for (j = 0; j < p; j++)
      mu[l + j*k] = 0.;
  }

  // TODO (josh 2021/09/05): to get rid of the errors when using NA values, have
  // to change how mu is calculted, make it such that each component gets
  // divided by the number of non-missing values in the cluster (for that
  // component).
  for (i = 0; i < n; i++)
  {
    ii = ic1[i];
    for (j = 0; j < p; j++)
      mu[ii + j*k] += x[i + j*n];
  }

  for (j = 0; j < p; j++)
  {
    for (l = 0; l < k; l++)
      mu[l + j*k] /= (double) n_k[l];

    /* for (i = 0; i < n; i++) */
    /* { */
    /*   ii = ic1[i]; */
    /*   da = x[i + j*n] - mu[ii + j*k]; */
    /*   sse[ii] += (da * da); */
    /* } */
  }

  return npass;
}

/*
 * Function to perform the optimal transfer step of the Hartigan-Wong algorithm.
 * Appended a k on the end to be able to test multiple version of H-W algorithm
 * in one testing program
 *
 * @param *x             - n*p length array containing col-major data
 * @param *mu            - k*p array of cluster centers
 * @param *sse           - k array of SSE for each cluster
 * @param *an1           - k work array for internal use an1[l] = n_k[l] / (n_k[l] - 1)
 * @param *an2           - k work array for internal use an2[l] = n_k[l] / (n_k[l] + 1)
 * @param *n_k           - k array containing number in each cluster
 * @param  n             - number of data points
 * @param  p             - dimension of data
 * @param  k             - number of clusters
 * @param *ic1           - n array of closest cluster indicator
 * @param *ic2           - n array of 2nd closest cluster indicator
 * @param *ncp           - k work array
 * @param *itran         - k work array
 * @param *indx          - 1 work array (only an array bc of pointer issue)
 * @param *d             - n work array
 * @param *live          - k work array (live set)
 * @param *kern_cross    - k array containing (sum i)(sum j)K(x_i, x_j) fe clust
 * @param *fo1           - n work array
 * @param *fo2           - n work array
 */
int optimal_transfer(double *x,
                     double *mu,
                     double *sse,
                     double *n_minus,
                     double *n_plus,
                     int    *n_k,
                     int     n,
                     int     p,
                     int     k,
                     int    *ic1,
                     int    *ic2,
                     int    *change,
                     int    *itran,
                     double *loss,
                     int    *live,
                     double *kern_cross,
                     double *kernel_matrix,
                     double *fo1,
                     double *fo2)
{
  /* c1 is the current closest cluster, c2 is the current second closest, and cl
   * is the second closest for an observation upon entering this stage */
  int c1, c2;

  int i, j, l;

  int flag = 0;

  double M, temp;

  double acc, fo_acc, self_kern;

  double big = 1.0E+10;

  for (l = 0; l < k; l++)
  {
    if (itran[l])
      live[l] = n + 1;
  }

  for (i = 0; i < n; i++)
  {
    self_kern = kernel_matrix[get_index(i, i, n)];

    c1 = ic1[i];
    c2 = ic2[i];

    /* skip this observation if it is the only member of its cluster */
    if (n_k[c1] == 1) continue;

    /* if cluster c1 has been updated in this stage, we need to recompute the
     * loss for this observation */
    if (change[c1])
    {
      acc = self_kern + 1. / (n_k[c1] * n_k[c1]) * kern_cross[c1];
      fo_acc = 0;

      for (j = 0; j < n; j++)
        if (ic1[j] == c1) fo_acc += kernel_matrix[get_index(i, j, n)];

      fo1[i] = fo_acc;
      acc -= 2. / (n_k[c1]) * fo_acc;

      loss[i] = acc * n_minus[c1];
    }

    /* this is the same thing as above, except with c2 */
    if (change[c2])
    {
      fo_acc = 0;

      for (j = 0; j < n; j++)
        if (ic1[j] == c2) fo_acc += kernel_matrix[get_index(i, j, n)];

      fo2[i] = fo_acc;
    }

    /* acc = self_kern + 1. / (n_k[c2] * n_k[c2]) * kern_cross[c2]; */
    /* acc -= 2. / (n_k[c2]) * fo2[i]; */
    /* M = acc * n_plus[c2]; */
    M = big;

    for (l = 0; l < k; l++)
    {
      if (l == c1) continue;

      /* check if this cluster has the minumum if either c1 or l is in the live
       * set */
      if (live[c1] > i || live[l] > i)
      {
        acc = self_kern + 1. / (n_k[l] * n_k[l]) * kern_cross[l];

        fo_acc = 0;
        for (j = 0; j < n; j++)
          if (ic1[j] == l) fo_acc += kernel_matrix[get_index(i, j, n)];

        acc -= 2. / (n_k[l]) * fo_acc;

        if (acc * n_plus[l] < M)
        {
          M = acc * n_plus[l];
          c2 = l;
          fo2[i] = fo_acc;
        }
      }
    }

    if (loss[i] <= M)
      ic2[i] = c2;
    else
    {
      flag++;
      live[c1] = n + i;
      live[c2] = n + i;

      itran[c1] = 1;
      itran[c2] = 1;

      change[c1] = i;
      change[c2] = i;


      kern_cross[c1] -= 2*fo1[i] + self_kern;
      kern_cross[c2] += 2*fo2[i] + 2*self_kern;

      n_k[c1]--;
      n_k[c2]++;

      if (n_k[c1] <= 1)
        n_minus[c1] = big;
      else
      {
        n_minus[c1] = (double) n_k[c1];
        n_minus[c1] /= (n_k[c1] - 1);
      }

      n_plus[c1] = (double) n_k[c1];
      n_plus[c1] /= (n_k[c1] + 1);

      if (n_k[c2] <= 1)
        n_minus[c2] = big;
      else
      {
        n_minus[c2] = (double) n_k[c2];
        n_minus[c2] /= (n_k[c2] - 1);
      }

      n_plus[c2] = (double) n_k[c2];
      n_plus[c2] /= (n_k[c2] + 1);

      ic1[i] = c2;
      ic2[i] = c1;

      temp = fo2[i];
      fo2[i] = fo1[i] - self_kern;
      fo1[i] = temp + self_kern;
    }

    /* if (indx[0] == n) */
    /* { */
    /*   return; */
    /* } */
  }

  for (l = 0; l < k; l++)
    live[l] -= n;

  return flag;
}

/*
 * Function to perform the quick transfer step of the Hartigan-Wong algorithm.
 * Appended a k on the end to be able to test multiple version of H-W algorithm
 * in one testing program
 *
 * @param *x          - n*p length array containing col-major data
 * @param *mu         - k*p array of cluster centers
 * @param *sse        - k array of SSE for each cluster
 * @param *an1        - k work array for internal use
 * @param *an2        - k work array for internal use
 * @param *n_k        - k array containing number in each cluster
 * @param  n          - number of data points
 * @param  p          - dimension of data
 * @param  k          - number of clusters
 * @param *ic1        - n array of closest cluster indicator
 * @param *ic2        - n array of 2nd closest cluster indicator
 * @param *ncp        - k work array
 * @param *itran      - k work array
 * @param *indx       - 1 work array (only an array bc of pointer issue)
 * @param *d          - n work array
 * @param *live       - k work array (live set)
 * @param *kern_cross - k array containing (sum i)(sum j)K(x_i, x_j) fe clust
 * @param *fo1        - n work array
 * @param *fo2        - n work array
 */
void quick_transfer(double *x,
                    double *mu,
                    double *sse,
                    double *n_minus,
                    double *n_plus,
                    int    *n_k,
                    int     n,
                    int     p,
                    int     k,
                    int    *ic1,
                    int    *ic2,
                    int    *change,
                    int    *itran,
                    double *loss,
                    int    *live,
                    double *kern_cross,
                    double *kernel_matrix,
                    double *fo1,
                    double *fo2)
{
  int c1, c2;
  int i, j;
  double big = 1.0E+10;
  double self_kern;

  int flag, nstep;

  double acc, fo_acc, temp;

  int *ls = (int *) S_alloc(k, sizeof(int));

  for (int l = 0; l < k; l++)
  {
    if (itran[l])
    {
      ls[l] = 1;
      itran[l] = 0;
    }
  }

  nstep = 0;

  double old_sse = big;
  double epsilon = .01;

  while (1)
  {
    flag = 0;

    if (nstep > 0)
    {
      old_sse = 0;
      for (i = 0; i < n; i++)
        old_sse += loss[i];
    }

    for (i = 0; i < n; i++)
    {
      c1 = ic1[i];
      c2 = ic2[i];
      nstep++;

      self_kern = kernel_matrix[get_index(i, i, n)];

      if (n_k[c1] == 1) continue;

      if (nstep <= change[c1])
      {
        acc = self_kern + 1. / (n_k[c1] * n_k[c1]) * kern_cross[c1];
        fo_acc = 0;

        for (j = 0; j < n; j++)
          if (ic1[j] == c1) fo_acc += kernel_matrix[get_index(i, j, n)];

        fo1[i] = fo_acc;
        acc -= 2. / (n_k[c1]) * fo_acc;

        loss[i] = acc * n_minus[c1];
      }

      if (nstep <= change[c2])
      {
        fo_acc = 0;

        for (j = 0; j < n; j++)
          if (ic1[j] == c2) fo_acc += kernel_matrix[get_index(i, j, n)];

        fo2[i] = fo_acc;
      }

      if (!ls[c1] && !ls[c2]) continue;

      acc = self_kern + 1. / (n_k[c2] * n_k[c2]) * kern_cross[c2];
      fo_acc = fo2[i];
      acc -= 2. / (n_k[c2]) * fo_acc;

      if (acc * n_plus[c2] < loss[i])
      {
        flag = 1;

        itran[c1] = 1;
        itran[c2] = 1;

        change[c1] = nstep + n;
        change[c2] = nstep + n;

        kern_cross[c1] -= 2*fo1[i] + self_kern;
        kern_cross[c2] += 2*fo2[i] + self_kern;

        n_k[c1]--;
        n_k[c2]++;

        if (n_k[c1] <= 1)
          n_minus[c1] = big;
        else
        {
          n_minus[c1] = (double) n_k[c1];
          n_minus[c1] /= (n_k[c1] - 1);
        }

        n_plus[c1] = (double) n_k[c1];
        n_plus[c1] /= (n_k[c1] + 1);

        if (n_k[c2] <= 1)
          n_minus[c2] = big;
        else
        {
          n_minus[c2] = (double) n_k[c2];
          n_minus[c2] /= (n_k[c2] - 1);
        }

        n_plus[c2] = (double) n_k[c2];
        n_plus[c2] /= (n_k[c2] + 1);

        ic1[i] = c2;
        ic2[i] = c1;

        temp = fo2[i];
        fo2[i] = fo1[i] - self_kern;
        fo1[i] = temp + self_kern;

      }
    }

    acc = 0;
    for (i = 0; i < n; i++)
      acc += loss[i];

    // Rprintf("%f\n", (old_sse - acc) / old_sse);
    // Rprintf("%f\n", old_sse);

    if (!flag || (old_sse - acc) / old_sse < epsilon) return;
  }
}

int get_index(int i, int j, int n)
{
  return (i + j*n);
}

int rand_dunif(int r)
{
  GetRNGstate();
  double u = unif_rand();
  PutRNGstate();

  return (trunc(u * r));
}

int rand_multinom(int n, double *probs)
{
  /* make the cdf */
  double *cdf = (double *) S_alloc(n + 1, sizeof(double));

  GetRNGstate();
  double u = unif_rand();
  PutRNGstate();

  for (int i = 1; i <= n; i++)
    cdf[i] = cdf[i - 1] + probs[i];

  for (int i = 1; i <= n; i++)
  {
    if (u <= cdf[i] && u > cdf[i - 1])
    {
      return i;
    }
  }
  /* this should never happen */
  return -1;
}

// TODO: these shouldn't be inthis file?
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

void center_kernel_matrix(double *kernel_matrix, int n)
{
  double *rowmeans = (double *) S_alloc(n, sizeof(double));
  double acc;
  double summean = 0;
  for (int i = 0; i < n; i++)
  {
    acc = 0;
    for(int j = 0; j < n; j++)
      acc += kernel_matrix[get_index(i, j, n)];
    rowmeans[i] = acc / n;
    summean += rowmeans[i];
  }
  summean /= n;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++)
      kernel_matrix[get_index(i, j, n)] += summean - rowmeans[i] - rowmeans[j];
  return;
}


int macqueen_step(double *x,
                  double *mu,
                  double *sse,
                  double *n_minus,
                  double *n_plus,
                  int    *n_k,
                  int     n,
                  int     p,
                  int     k,
                  int    *ic1,
                  double *loss,
                  double *kern_cross,
                  double *kernel_matrix,
                  double *fo1)
{

  /* c1 is the current closest center, c2 is the candidate cluster */
  int c1, c2;
  int i, j, l;
  int n_transfer = 0;
  double acc;

  double self_kern;

  double M, temp, fo_acc, fo_low;
  double big = 1.0E+10;


  for (i = 0; i < n; i++)
  {
    self_kern = kernel_matrix[get_index(i, i, n)];

    c1 = ic1[i];

    /* skip if the observation is the only member of its cluster */
    if (n_k[c1] == 1) continue;

    acc = self_kern + 1. / (n_k[c1] * n_k[c1]) * kern_cross[c1];
    fo_acc = 0;

    for (j = 0; j < n; j++)
      if (ic1[j] == c1) fo_acc += kernel_matrix[get_index(i, j, n)];

    fo1[i] = fo_acc;
    acc -= 2. / (n_k[c1]) * fo_acc;

    loss[i] = acc * n_minus[c1];

    M = big;
    c2 = -1;

    for (l = 0; l < k; l++)
    {
      if (l == c1) continue;
      acc = self_kern + 1. / (n_k[l] * n_k[l]) * kern_cross[l];

      fo_acc = 0;
      for (j = 0; j < n; j++)
        if (ic1[j] == l) fo_acc += kernel_matrix[get_index(i, j, n)];

      acc -= 2. / (n_k[l]) * fo_acc;

      if (acc < M)
      {
        M = acc;
        c2 = l;
        fo_low = fo_acc;
      }
    }
    if (M < loss[i] / n_minus[c1])
    {
      n_transfer++;

      kern_cross[c1] -= 2*fo1[i] + self_kern;
      kern_cross[c2] += 2*fo_low + 2*self_kern;

      n_k[c1]--;
      n_k[c2]++;

      if (n_k[c1] <= 1)
        n_minus[c1] = big;
      else
      {
        n_minus[c1] = (double) n_k[c1];
        n_minus[c1] /= (n_k[c1] - 1);
      }

      n_plus[c1] = (double) n_k[c1];
      n_plus[c1] /= (n_k[c1] + 1);

      if (n_k[c2] <= 1)
        n_minus[c2] = big;
      else
      {
        n_minus[c2] = (double) n_k[c2];
        n_minus[c2] /= (n_k[c2] - 1);
      }

      n_plus[c2] = (double) n_k[c2];
      n_plus[c2] /= (n_k[c2] + 1);

      ic1[i] = c2;
      fo1[i] = fo_low + self_kern;
    }
  }

  return n_transfer;
}
int lloyd_step(double *x,
               double *mu,
               double *sse,
               double *n_minus,
               double *n_plus,
               int    *n_k,
               int     n,
               int     p,
               int     k,
               int    *ic1,
               double *loss,
               double *kern_cross,
               double *kernel_matrix,
               double *fo1)
{

  /* c1 is the current closest center, c2 is the candidate cluster */
  int c1, c2;
  int i, j, l;
  int n_transfer = 0;
  double acc;
  int *new_clusters = (int *) S_alloc(n, sizeof(int));
  memcpy(new_clusters, ic1, n*sizeof(int));

  double self_kern;

  double M, temp, fo_acc, fo_low;
  double big = 1.0E+10;


  for (i = 0; i < n; i++)
  {
    self_kern = kernel_matrix[get_index(i, i, n)];

    c1 = ic1[i];

    /* skip if the observation is the only member of its cluster */
    if (n_k[c1] == 1) continue;

    acc = self_kern + 1. / (n_k[c1] * n_k[c1]) * kern_cross[c1];
    fo_acc = 0;

    for (j = 0; j < n; j++)
      if (ic1[j] == c1) fo_acc += kernel_matrix[get_index(i, j, n)];

    fo1[i] = fo_acc;
    acc -= 2. / (n_k[c1]) * fo_acc;

    loss[i] = acc * n_minus[c1];

    M = big;
    c2 = -1;

    for (l = 0; l < k; l++)
    {
      if (l == c1) continue;
      acc = self_kern + 1. / (n_k[l] * n_k[l]) * kern_cross[l];

      fo_acc = 0;
      for (j = 0; j < n; j++)
        if (ic1[j] == l) fo_acc += kernel_matrix[get_index(i, j, n)];

      acc -= 2. / (n_k[l]) * fo_acc;

      if (acc < M)
      {
        M = acc;
        c2 = l;
        fo_low = fo_acc;
      }
    }
    if (M < loss[i] / n_minus[c1])
    {
      n_transfer++;
      new_clusters[i] = c2;
    }
  }

  memcpy(ic1, new_clusters, n*sizeof(int));
  return n_transfer;
}

void init_centers(int n, int k, int *ic1, int *ic2, int *n_k, double *kernel_matrix)
{
  int i, j;
  double self_kern = kernel_matrix[get_index(i, i, n)];
  int first_center = rand_dunif(n);
  int new_center, old_clust;
  double *closest_dist = (double *) S_alloc(n, sizeof(double));
  double *new_dist = (double *) S_alloc(n, sizeof(double));

  double aux;

  double total_prob;
  double *prob_vec = (double *) S_alloc(n, sizeof(double));
  double *cum_sum = (double *) S_alloc(n, sizeof(double));

  for (i = 0; i < n; i++)
    ic1[i] = 0;

  /* now compute the distance from each point to the cluster with the single
   * point */
  for (i = 0; i < n; i++)
    closest_dist[i] = 2*self_kern - 2*kernel_matrix[get_index(i, first_center, n)];

  total_prob = 0;
  for (i = 0; i < n; i++)
  {
    total_prob += closest_dist[i];
    cum_sum[i] = total_prob;
  }

  for (i = 0; i < n; i++)
    prob_vec[i] = closest_dist[i] / total_prob;

  for (j = 1; j < k; j++)
  {
    GetRNGstate();
    aux = unif_rand();
    PutRNGstate();

    for (i = 0; i < n; i++) {
      if (cum_sum[i] - aux >= 0) break;
    }
    /* i is now the chosen point for the new cluster */
    new_center = i;
    ic1[new_center] = j;

    for (i = 0; i < n; i++)
    {
      new_dist[i] = 2*self_kern - 2*kernel_matrix[get_index(i, new_center, n)];
      if (new_dist[i] < closest_dist[i])
      {
        closest_dist[i] = new_dist[i];
        old_clust = ic1[i];
        ic1[i] = j;
        ic2[i] = old_clust;
      }
    }
    total_prob = 0;
    for (i = 0; i < n; i++)
    {
      total_prob += closest_dist[i];
      cum_sum[i] = total_prob;
    }
  }

  for (j = 0; j < k; j++)
    for (i = 0; i < n; i++)
      if (ic1[i] == j)
        n_k[j]++;

  return;
}
