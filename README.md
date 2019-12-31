
# Fast Mutual Information based independence Test (fastmit)

### Introduction
The fundamental problem for data mining and statistical analysis is:

- Whether two random variables are dependent?

**fastmit** package provides solutions for this issue. It implements the kNN method described by Kraskov, et. al (2004) to estimate the empirical mutual information and furthermore uses permutation test to detect whether two random variables are independent. The core functions in **fastmit** package are **mi** and **mi.test**.

These functions based on mutual information have two main advantages:

- It's applicable to complex data in metric spaces.

- It is faster than other statistics (e.g., distance covariance and ball covariance), which makes it advantageous in large sample situations.

## Installation
#### CRAN version
You can install the released version of fastmit from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("fastmit")
```

#### Github version
You can install the released version of fastmit from GitHub with:

``` r
library(devtools)
install_github("Mamba413/fastmit")
```
*Windows* user will need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first.    



## Example
- **mi** function
``` r
## simulate data
set.seed(1)
x <- rnorm(100)
y <- x + rnorm(100)
## estimate the empirical mutual information
mi(x, y)
```
In this example, the result is:
```
# [1] 0.3320034
```

- **mi.test** function

``` r
## simulate data
set.seed(1)
error <- runif(50, min = -0.3, max = 0.3)
x <- runif(50, 0, 4*pi)
y <- cos(x) + error

## perform independence test via mutual information
mi.test(x, y)
```
In this example, the result is:
```
	Mutual Information test of independence

data:  x and y
number of observations = 50
replicates = 99
p-value = 0.01
alternative hypothesis: random variables are dependent
sample estimates:
       MI 
0.6953105 
```

If you find any bugs, or if you experience any crashes, please report to us. Also, if you have any questions, feel free to ask.

### License
GPL-3

### Reference
- Alexander Kraskov; Harald Stögbauer; Peter Grassberger
Phys. [Estimating Mutual Information](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.69.066138) Rev. E 69, 066138
