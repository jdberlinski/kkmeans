/**
 * @file kkmeans.c
 * @auth Josh Berlinski
 * @date (last edit) 2020/05/18: moved all of the kernel stuff to a separate file, removed the
 * read_tabular function (unneeded in here anyway)
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

/* perhaps I should add a header file by now... */
int rand_dunif(int r);
int get_index(int i, int j, int n);
void optrak(double x[], double mu[], double sse[], double an1[], double an2[], int n_k[], int n, int p,
    int k, int ic1[], int ic2[], int ncp[], int itran[], int indx[], double d[], int live[],
    double kern_cross[], double kernel_matrix[], double fo1[], double fo2[], double h);
void qtrank(double x[], double mu[], double sse[], double an1[], double an2[], int n_k[], int n, int p,
    int k, int ic1[], int ic2[], int ncp[], int itran[], int indx[], double d[], int live[],
    double kern_cross[], double kernel_matrix[], double fo1[], double fo2[], double h);
/*
 * The main clustering algorithm, with the goal of dividing n observations into
 * k clusters such that the WSS (in feature space) is minimized
 *
 * @param x[]      - n*p length array containing col-major data
 * @param n        - number of data points
 * @param p        - dimension of data
 * @param k        - number of clusters
 * @param h        - tuning parameter for kernel function
 * @param iter_max - the maximum number of passes allowed
 * @param kernel   - the kernel function to use
 * @param mu[]     - k*p array of cluster centers
 * @param sse[]    - k array of SSE for each cluster
 * @param ic1[]    - n array of closest cluster indicator
 */
void kcluster(double x[], int n, int p, int k, double h, int iter_max,
    double (*kernel)(int, int, double[], int, int, double),
    double mu[], double sse[], int ic1[])
{
  // I'm not sure if this is really necessary, but it's the safest option, at
  // least. R will allocate the memory (and initialize to 0, then clear it after
  double* an1        = (double *) S_alloc(k, sizeof(double));
  double* an2        = (double *) S_alloc(k, sizeof(double));
  double* d          = (double *) S_alloc(n, sizeof(double));
  double* fo1        = (double *) S_alloc(n, sizeof(double));
  double* fo2        = (double *) S_alloc(n, sizeof(double));
  double* kern_cross = (double *) S_alloc(k, sizeof(double));

  int* ic2   = (int *) S_alloc(n, sizeof(int));
  int* ncp   = (int *) S_alloc(k, sizeof(int));
  int* itran = (int *) S_alloc(k, sizeof(int));
  int* live  = (int *) S_alloc(k, sizeof(int));
  int* n_k   = (int *) S_alloc(k, sizeof(int));

  double* kernel_matrix = (double *) S_alloc( n*n , sizeof(double));
  double big = 1.0E+10;

  int i, j, l, ij, ii, iw, iu;
  double da, aa;

  int indx[1] = { 0 };
  int current_clust = 0;

  /* original cluster assignment here */

  /* randomize the clusters */
  for( int j = 0; j < n; j++ )
  {
    current_clust = rand_dunif(k);

    n_k[current_clust]++;

    ic1[j] = current_clust;
    ic2[j] = (current_clust + 1) % k;
  }


  for( i = 0; i < n; i++ )
    for( j = 0; j <= i; j++ )
    {
      kernel_matrix[get_index(i, j, n)] = kernel(i, j, x, n, p, h);
      kernel_matrix[get_index(j, i, n)] = kernel(i, j, x, n, p, h);
    }

  /* choose the first center at random among the data points */
  /* D(x) = ||x - c||^2
   *      = (x - c)'(x - c)
   *      = (x' - c')(x - c)
   *      = x'x - x'c - c'x  + c'c
   *      = k(x, x) - 2k(x, c) + k(c, c)
   */


  /* initialize the last term */
  for( l = 0; l < k; l++ )
  {
    kern_cross[l] = 0.;

    for( iw = 0; iw < n; iw++ )
      if( ic1[iw] == l )
        for( iu = 0; iu < n; iu++ )
          if( ic1[iu] == l )
            kern_cross[l] += kernel_matrix[get_index(iw, iu, n)];
  } // end for l

  /* update cluster centers to be the average of points contained within them */
  for( l = 0; l < k; l++ )
  {
    n_k[l] = 0;
    for( j = 0; j < p; j++)
      mu[l + j*k] = 0.;
  } // end for l

  for( i = 0; i < n; i++ )
  {
    l = ic1[i];
    n_k[l]++;
    for( j = 0; j < p; j++ )
      mu[l + j*k] += x[i + j*n];
  } // end for i

  for( l = 0; l < k; l++ )
    if( n_k[l] == 0 )
      error("Error: Cluster %d has 0 observations. Exiting...\n Possible bad starting seed.\n", l);

  for( l = 0; l < k; l++ )
  {
    aa = (double) n_k[l];

    for( j = 0; j < p; j++ )
      mu[l + j*k] /= aa;

    /* initialize an1, an2, itran, ncp
     * an1[l] = n_k[l] / (n_k[l] - 1)
     * an2[l] = n_k[l] / (n_k[l] + 1)
     * itran[l] = 1 if cluster l is updated in the qtran step
     *          = otherwise*/

    an2[l] = aa / (aa + 1.);
    if( aa > 1. )
      an1[l] = aa / (aa - 1.);
    else
      an1[l] = big;

    itran[l] = 1;
    ncp[l] = -1;

  } // end for l

  indx[0] = 0;

  for( ij = 0; ; ij++ )
  {
    /* in this stage, there is only one pass through the data. each point is
     * re-allocated, if necessary, to the cluster that will induce the maximum
     * reduction in within-cluster sum of squares */

    /* printf("optra\n"); */
    optrak(x, mu, sse, an1, an2, n_k, n, p, k, ic1, ic2, ncp, itran, indx, d, live,
        kern_cross, kernel_matrix, fo1, fo2, h);

    /* stop if no transfer took place in the last n optimal transfer steps */
    if( indx[0] == n ) // goto 150? (later)
      break;

    /* printf("qtran\n"); */
    qtrank(x, mu, sse, an1, an2, n_k, n, p, k, ic1, ic2, ncp, itran, indx, d, live,
        kern_cross, kernel_matrix, fo1, fo2, h);

    /* if there are only two clusters, there is no need to re-enter the optran
     * stage */
    if( k == 2 ) //goto 150??? (later)
      break;
    /* ncp has to be set to 0 before entering optra */
    for( l = 0; l < k; l++ )
      ncp[l] = 0;

    if( ij >= iter_max)
      break;
  } // end for ij

  /* ifault error checking goes here */
  /* compute the within-cluster ss for each clustr */
  for( l = 0; l < k; l++ )
  {
    sse[l] = 0.;
    for( j = 0; j < p; j++ )
      mu[l + j*k] = 0.;
  } // end for l

  for( i = 0; i < n; i++ )
  {
    ii = ic1[i];
    for( j = 0; j < p; j++ )
      mu[ii + j*k] += x[i + j*n];
  } // end for i

  for( j = 0; j < p; j++ )
  {
    for( l = 0; l < k; l++ )
      mu[l + j*k] /= (double) n_k[l];

    for( i = 0; i < n; i++ )
    {
      ii = ic1[i];
      da = x[i + j*n] - mu[ii + j*k];
      sse[ii] += (da * da);
    } // end for i
  } // end for j
  return;
}

/*
 * The kernel function, just a standard Gaussian kernel for now.
 * TODO: make this into a function pointer?
 *
 * @param norm the value to compute the kernel at (ie ||x_i - x_j||^2)
 * @param h tuning parameter
 * @return K(||x_i - x_j||^2)
 */
//double kernel(double norm, double h)
//{
 // return( exp(-0.5 * norm / h) );
//}



/*
 * Function to perform the optimal transfer step of the Hartigan-Wong algorithm.
 * Appended a k on the end to be able to test multiple version of H-W algorithm
 * in one testing program
 *
 * @param x[]           - n*p length array containing col-major data
 * @param mu[]          - k*p array of cluster centers
 * @param sse[]         - k array of SSE for each cluster
 * @param an1[]         - k work array for internal use an1[l] = n_k[l] / (n_k[l] - 1)
 * @param an2[]         - k work array for internal use an2[l] = n_k[l] / (n_k[l] + 1)
 * @param n_k[]         - k array containing number in each cluster
 * @param n             - number of data points
 * @param p             - dimension of data
 * @param k             - number of clusters
 * @param ic1[]         - n array of closest cluster indicator
 * @param ic2[]         - n array of 2nd closest cluster indicator
 * @param ncp[]         - k work array
 * @param itran[]       - k work array
 * @param indx[]        - 1 work array (only an array bc of pointer issue)
 * @param d[]           - n work array
 * @param live[]        - k work array (live set)
 * @param kern_cross[]  - k array containing (sum i)(sum j)K(x_i, x_j) fe clust
 * @param fo1[]         - n work array
 * @param fo2[]         - n work array
 */
void optrak(double x[], double mu[], double sse[], double an1[], double an2[], int n_k[], int n, int p,
    int k, int ic1[], int ic2[], int ncp[], int itran[], int indx[], double d[], int live[],
    double kern_cross[], double kernel_matrix[], double fo1[], double fo2[], double h)
{

  int l, i, l1, l2, ll, j, al1, al2, iu;
  double da, dc, de, M, alt, alw;
  double fo_kern, self_kern, temp;
  double big = 1.0E+10;
  /* If cluster l is updated in the last qtran stage, it belongs to the live set
   * throughout this stage, otherwise at each step  it is not in the live set if
   * it has not been updated in the last M optra steps*/
  for( l = 0; l < k; l++ )
    if( itran[l] == 1 ) live[l] = n + 1;

  for( i = 0; i < n; i++ )
  {
    self_kern = kernel_matrix[get_index(i, i, n)];
    indx[0]++;
    l1 = ic1[i];
    l2 = ic2[i];
    ll = l2;

    /* if point i is the only member of cluster l1, no transfer */
    if( n_k[l1] > 1 )
    {
      /* if l1 has not yet been updated in this stage, no need to recompute d[i] */
      if( ncp[l1] != 0 )
      {
        de = self_kern + 1. / (n_k[l1] * n_k[l1]) * kern_cross[l1];
        fo_kern = 0.;
        for( iu = 0; iu < n; iu++ )
          if( ic1[iu] == l1 )
            fo_kern += kernel_matrix[get_index(i, iu, n)];
        fo1[i] = fo_kern;
        de -= 2. / (n_k[l1]) * fo_kern;
        d[i] = de * an1[l1];
      } // end if

      /* THis might not be necessary, but check if l2 has been updated at this
       * stage, and recompute (might be a big brain move though) */
      if( ncp[l2] != 0 )
      {
        fo_kern = 0.;
        for( iu = 0; iu < n; iu++ )
          if( ic1[iu] == l2 )
            fo_kern += kernel_matrix[get_index(i, iu, n)];
        fo2[i] = fo_kern;
      } // end if

      /* find the cluster with minimum M*/

      da = self_kern + 1. / (n_k[l2] * n_k[l2]) * kern_cross[l2];
      fo_kern = fo2[i];
      da -= 2. / (n_k[l2]) * fo_kern;
      M = da * an2[l2];


      for( l = 0; l < k; l++ )
      {
        /* if live[l1] <= i, then l1 is not in the live set.
         * f this is true, we only need to consider clusters that are in the
         * live set for possible transfers of point i. Otherwise we need to
         * consider all possible clusters */
        if( (i < live[l1] || i < live[l2]) && (l != l1 && l != ll) )
        {
          dc = self_kern + 1. / (n_k[l] * n_k[l]) * kern_cross[l];
          fo_kern = 0.;
          for( iu = 0; iu < n; iu++ )
            if( ic1[iu] == l )
              fo_kern += kernel_matrix[get_index(i, iu, n)];
          dc -= 2. / (n_k[l]) * fo_kern;

          if( dc*an2[l] < M )
          {
            M = dc * an2[l];
            l2 = l;
            fo2[i] = fo_kern;
          } // end if
        } // end if
      } // end for l

      /* if no transfer is necesary, l2 is the new ic2 */
      if( d[i] <= M )
      {
        ic2[i] = l2;
      } else
      {
        indx[0] = 0;
        live[l1] = n + i;
        live[l2] = n + i;
        ncp[l1] = i;
        ncp[l2] = i;
        al1 = n_k[l1];
        alw = al1 - 1.; // alw = n-1
        al2 = n_k[l2];
        alt = al2 + 1.; // alt = n / n+1

        /* the c computation for updating means is the same */
        for( j = 0; j < p; j++ )
        {
          mu[l1 + j*k] = (mu[l1 + j*k] * al1 - x[i + j*n]) / alw;
          mu[l2 + j*k] = (mu[l2 + j*k] * al2 + x[i + j*n]) / alt;
        } // end for j

        n_k[l1]--;
        n_k[l2]++;
        al1 = n_k[l1];
        alw = al1 - 1.;
        al2 = n_k[l2];
        alt = al2 + 1.;
        kern_cross[l1] -= 2*fo1[i] + self_kern;
        kern_cross[l2] += 2*fo2[i] + 2*self_kern;
        an2[l1] = alw / al1;

        if( alw > 1. )
          an1[l1] = alw / (alw - 1.);
        else
          an1[l1] = big;

        an1[l2] = alt / al2;
        an2[l2] = alt / (alt + 1.);

        ic1[i] = l2;
        ic2[i] = l1;

        temp = fo2[i];
        fo2[i] = fo1[i] - self_kern;
        fo1[i] = temp + self_kern;
      } // end else
    } // end if

    if( indx[0] == n )
    {
      return;
    } // end if
  } // end for i

  for( l = 0; l < k; l++ )
  {
    itran[l] = 0;
    live[l] -= n;
  } // end for l

  return;
}

/*
 * Function to perform the quick transfer step of the Hartigan-Wong algorithm.
 * Appended a k on the end to be able to test multiple version of H-W algorithm
 * in one testing program
 *
 * @param x[]           - n*p length array containing col-major data
 * @param mu[]          - k*p array of cluster centers
 * @param sse[]         - k array of SSE for each cluster
 * @param an1[]         - k work array for internal use
 * @param an2[]         - k work array for internal use
 * @param n_k[]         - k array containing number in each cluster
 * @param n             - number of data points
 * @param p             - dimension of data
 * @param k             - number of clusters
 * @param ic1[]         - n array of closest cluster indicator
 * @param ic2[]         - n array of 2nd closest cluster indicator
 * @param ncp[]         - k work array
 * @param itran[]       - k work array
 * @param indx[]        - 1 work array (only an array bc of pointer issue)
 * @param d[]           - n work array
 * @param live[]        - k work array (live set)
 * @param kern_cross[]  - k array containing (sum i)(sum j)K(x_i, x_j) fe clust
 * @param fo1[]         - n work array
 * @param fo2[]         - n work array
 */
void qtrank(double x[], double mu[], double sse[], double an1[], double an2[], int n_k[], int n, int p,
    int k, int ic1[], int ic2[], int ncp[], int itran[], int indx[], double d[], int live[],
    double kern_cross[], double kernel_matrix[], double fo1[], double fo2[], double h)
{
  /* in the optra stage, ncp indicates the step at which cluster l is last
   * udpated. In the qtran stage, ncp is equal to the step at which cluster l is
   * last updated plus n */
  int icoun = 0, istep = 0;
  int i, l1, l2, j, iu;
  double da, dd, alt, alw, al1, al2;
  double fo_kern, self_kern, temp;
  double big = 1.0E+10;

  while( 1 )
  {

    for( i  = 0; i < n; i++ )
    {
      icoun++;
      istep++;


      l1 = ic1[i];
      l2 = ic2[i];

      self_kern = kernel_matrix[get_index(i, i, n)];

      /* if point i is the only member of cluster l1, no transfer */
      if( n_k[l1] > 1 )
      {
        /* if ncp(l1) < istep, no need to recompute the distance from point i to
         * cluster l1. note that if cluster l1 is the last update exactly n steps
         * ago, we still need to compute the distance from point i to cluster l1 */
        if( istep <= ncp[l1] )
        {
          da = self_kern + 1. / (n_k[l1] * n_k[l1]) * kern_cross[l1];
          fo_kern = 0.;
          for( iu = 0; iu < n; iu++ )
            if( ic1[iu] == l1 )
              fo_kern += kernel_matrix[get_index(i, iu, n)];
          fo1[i] = fo_kern;
          da -= 2. / (n_k[l1]) * fo_kern;
          d[i] = da * an1[l1];
        } // end if

        /* same thing as optra, this might be dumb */
        if( istep <= ncp[l2] )
        {
          fo_kern = 0.;
          for( iu = 0; iu < n; iu++ )
            if( ic1[iu] == l2 )
              fo_kern += kernel_matrix[get_index(i, iu, n)];
          fo2[i] = fo_kern;
        } // end if

        /* if ncp(l1) <= istep and ncp(l2) <= istep, there will be no transfer of
         * point i at this step */
        if( istep < ncp[l1] || istep < ncp[l2] )
        {
          dd = self_kern + 1. / (n_k[l2] * n_k[l2]) * kern_cross[l2];
          fo_kern = fo2[i];
          dd -= 2. / (n_k[l2]) * fo_kern;

          /* update cluster centers, ncp, nc, itran, an1, and an2, for cluster l1
           * and l2. also update ic1[i] and ic2[i].
           * note that if any updating occurs in this stage, indx is set bavk to 0 */
          if( dd*an2[l2] < d[i] )
          {
            icoun = 0;
            indx[0] = 0;
            itran[l1] = 1;
            itran[l2] = 1;
            ncp[l1] = istep + n;
            ncp[l2] = istep + n;

            al1 = n_k[l1];
            alw = al1 - 1.;
            al2 = n_k[l2];
            alt = al2 + 1.;
            for( j = 0; j < p; j++ )
            {
              mu[l1 + j*k] = ( mu[l1 + j*k] * al1 - x[i + j*n] ) / alw;
              mu[l2 + j*k] = ( mu[l2 + j*k] * al2 + x[i + j*n] ) / alt;
            } // end for j

            n_k[l1]--;
            n_k[l2]++;
            kern_cross[l1] -= 2*fo1[i] + self_kern;
            kern_cross[l2] += 2*fo2[i] + 2*self_kern;
            al1 = n_k[l1];
            alw = al1 - 1.;
            al2 = n_k[l2];
            alt = al2 + 1.;
            an2[l1] = alw / al1;
            if( alw > 1. )
              an1[l1] = alw / (alw - 1.);
            else
              an1[l1] = big;
            an1[l2] = alt / al2;
            an2[l2] = alt / (alt + 1.);
            ic1[i] = l2;
            ic2[i] = l1;
            temp = fo2[i];
            fo2[i] = fo1[i] - self_kern;
            fo1[i] = temp + self_kern;
          } //end if
        } // end if
      } // end if

      /* if no re-allocation took place in the last n steps, return */

      if( icoun == n )
      {
        return;
      }
    } // end for i
  } // end while( 1 )
}

int get_index(int i, int j, int n)
{
  return( i + j*n );
  /* return( i * (i+1)/2 +j ); */
}

int rand_dunif(int r)
{
  GetRNGstate();
  double u = unif_rand();
  PutRNGstate();

  return( trunc(u * r) );
}
