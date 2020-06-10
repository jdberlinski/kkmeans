library(kkmeans)
library(microbenchmark)

dpath <- "data/clustering_benchmarks/graves/"

ring <- read.table(paste0(dpath, "ring.data.gz"), header = FALSE)
line <- read.table(paste0(dpath, "line.data.gz"), header = FALSE)
ring_outliers <- read.table(paste0(dpath, "ring_outliers.data.gz"), header = FALSE)
zigzag <- read.table(paste0(dpath, "zigzag.data.gz"), header = FALSE)
para <- read.table(paste0(dpath, "parabolic.data.gz"), header = FALSE)

res_ring <- kkmeans(as.matrix(ring), 2, 'p', 2)
plot(ring, col = res_ring$cluster)

microbenchmark(
  kkmeans = kkmeans(as.matrix(ring), 2, 'p', 2),
  kernlab = kernlab::kkmeans(as.matrix(ring), 2, "polydot",
                             kpar = list(degree = 2, offset = 1, scale = 1))
)
# Unit: milliseconds
#     expr     min      lq    mean   median       uq      max neval
#  kkmeans 21.1545 23.4422  27.334  28.3861  30.5361  37.4312   100
#  kernlab 45.8327 86.4229 116.756 104.8772 139.0280 279.2102   100

res_zig <- kkmeans(as.matrix(zigzag), 3, 'g', .1)
plot(zigzag, col = res_zig$cluster)

res <- kernlab::kkmeans(as.matrix(zigzag), 3, "rbfdot", kpar = list(sigma = 10))
plot(zigzag, col = res@.Data)

microbenchmark(
  kkmeans = kkmeans(as.matrix(zigzag), 3, 'g', .1),
  kernlab = kernlab::kkmeans(as.matrix(zigzag), 3, "rbfdot", kpar = list(sigma = 10))
)
# Unit: milliseconds
#     expr      min      lq     mean   median       uq      max neval
#  kkmeans  2.20142 198.734 183.7727 200.0336 201.4637 209.7711   100
#  kernlab 16.70340  20.883  27.4017  25.3138  30.2379  84.7467   100

res_para <- kkmeans(as.matrix(para), 2, 'g', .7)
plot(para, col = res_para$cluster)
res_para <- kernlab::kkmeans(as.matrix(para), 2, "anovadot", kpar = list(sigma = 10, degree = 2))
plot(para, col = res_para@.Data)

