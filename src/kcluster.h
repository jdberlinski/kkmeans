void optrak(double *x,
            double *mu,
            double *sse,
            double *an1,
            double *an2,
            int    *n_k,
            int     n,
            int     p,
            int     k,
            int    *ic1,
            int    *ic2,
            int    *ncp,
            int    *itran,
            int    *indx,
            double *d,
            int    *live,
            double *kern_cross,
            double *kernel_matrix,
            double *fo1,
            double *fo2,
            double  h);

void qtrank(double *x,
            double *mu,
            double *sse,
            double *an1,
            double *an2,
            int    *n_k,
            int     n,
            int     p,
            int     k,
            int    *ic1,
            int    *ic2,
            int    *ncp,
            int    *itran,
            int    *indx,
            double *d,
            int    *live,
            double *kern_cross,
            double *kernel_matrix,
            double *fo1,
            double *fo2,
            double  h);

int get_index(int i, int j, int n);

int rand_dunif(int r);
int rand_multinom(int n, double *probs);

