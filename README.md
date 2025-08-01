# kkmeans

## Overview
The function `kkmeans::kkmeans(data, k, kern, params, iter_max)` performs kernel k-means using an algorithm
similar to [Hartigan and Wong (1979)][1]. See [Berlinski and Maitra (2025)][2]
for implementation details.

## Installation
In the extracted directory:
```{R}
devtools::install_github("jdberlinski/kkmeans")
```

## Usage
```{R}
library(kkmeans)

dat = iris[, 1:4]
res = kkmeans(as.matrix(dat), k = 3, kern = "gaussian", params = 1)
res
```

[1]: https://www.jstor.org/stable/2346830?seq=1
[2]: https://doi.org/10.1002/sam.70032
