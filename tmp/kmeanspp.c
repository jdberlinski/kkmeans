  /* choose the first center at random among the data points */
  /* D(x) = ||x - c||^2
   *      = (x - c)'(x - c)
   *      = (x' - c')(x - c)
   *      = x'x - x'c - c'x  + c'c
   *      = k(x, x) - 2k(x, c) + k(c, c)
   */

/*
 * This works in the sense that it doesn't give a segfault.
 * It doesn't work in the sense that it never gives the correct results.
 */

/* new cluster assignment using kmeans++ */

int init_clust = rand_dunif(n);
int* init_centers = (int *) S_alloc(k, sizeof(int));
init_centers[0] = init_clust;
double* dist_kpp = (double *) S_alloc(n, sizeof(double));
double* dist_kpp2 = (double *) S_alloc(n, sizeof(double));
int current_center = 0;

double t_dist = big;

for (i = 0; i < n; i++)
{
  /* ic1 is initially 0 */
  current_clust = ic1[i];
  current_center = init_centers[current_clust];
  dist_kpp[i] = kernel_matrix[get_index(i, i, n)] -
    2 * kernel_matrix[get_index(i, current_center, n)] +
    kernel_matrix[get_index(current_center, current_center, n)];

  dist_kpp[i] *= dist_kpp[i];

  dist_kpp2[i] = big;
  n_k[0]++;
}

for (j = 1; j < k; j++)
{
  /* TODO: either need to figure out how to get rmultinom to work or write a
   * faster in-house function */
  /* GetRNGstate(); */
  /* rmultinom(1, dist_kpp, n, &init_centers[j]); */
  /* PutRNGstate(); */
  init_centers[j] = rand_multinom(n, dist_kpp);

  current_clust = init_centers[j];

  for (i = 0; i < n; i++)
  {
    // t_dist is the distance from the point to its current center
    // but I already have that, so I need to compute the distance from each
    // point to the newly chosen center.
    // The distance from the point to its current center is going to be stored
    // in the dist_kpp matrix
    t_dist = kernel_matrix[get_index(i, i, n)] -
      2 * kernel_matrix[get_index(i, current_center, n)] +
      kernel_matrix[get_index(current_clust, current_clust, n)];

    t_dist *= t_dist;

    if (t_dist < dist_kpp[i])
    {
      ic2[i] = ic1[i];
      ic1[i] = current_clust;
      dist_kpp2[i] = dist_kpp[i];
      dist_kpp[i] = t_dist;
    }
    else if (t_dist < dist_kpp2[i])
    {
      ic2[i] = current_clust;
      dist_kpp2[i] = t_dist;
    }

  }
}

