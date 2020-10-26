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
                     double *fo2);
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
                    double *fo2);

int get_index(int i, int j, int n);

int rand_dunif(int r);
int rand_multinom(int n, double *probs);
void get_kernel_matrix(double *x,
                       int     n,
                       int     p,
                       double  h,
                       double  (*kernel)(int, int, double*, int, int, double),
                       double *kernel_matrix);
