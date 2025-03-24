library(MASS)
library(Matrix)
# define covariance matrix
orthog_mat <- function(n) {
  # orthogonal matrix
  # exponent <- mvrnorm(10, rep(0, 10), diag(10))
  # define a skew symmetric matrix
  A <- matrix(0, n, n)
  A[lower.tri(A)] <- rnorm(n * (n - 1) / 2)
  A <- A - t(A)
  return(matrix(expm(A), n, n))
}


create_cov <- function(n, alpha) {
  cov <- 0.5 * diag(10 * n)
  row_id_ <- 1
  col_id_ <- 11
  for (i in 1:(n - 1)) {
    row_id <- 10 * i
    col_id_ <- row_id_ + 10
    col_id <- row_id + 10
    for (j in (i + 1):n) {
      col_id <- 10 * j
      cov[row_id_:row_id, col_id_:col_id] <- alpha * orthog_mat(10)
      col_id_ <- col_id + 1
    }
    row_id_ <- row_id + 1
  }
  return(cov + t(cov))
}

create_cov_evals <- function(n, scalar) {
  cov <- 0.5 * diag(10 * n)
  row_id_ <- 1
  col_id_ <- 11
  for (i in 1:(n - 1)) {
    row_id <- 10 * i
    col_id_ <- row_id_ + 10
    col_id <- row_id + 10
    for (j in (i + 1):n) {
      col_id <- 10 * j
      mat <- orthog_mat(10) %*% diag(scalar * rep(1, 10)) %*% t(orthog_mat(10))
      cov[row_id_:row_id, col_id_:col_id] <- mat
      col_id_ <- col_id + 1
    }
    row_id_ <- row_id + 1
  }
  return(cov + t(cov))
}

create_cov_evals <- function(n, scalar) {
  cov <- 0.5 * diag(10 * n)
  row_id_ <- 1
  col_id_ <- 11
  row_id <- 10
  for (j in 2:n) {
    col_id <- 10 * j
    mat <- orthog_mat(10) %*% diag(scalar * rep(1, 10)) %*% t(orthog_mat(10))
    cov[row_id_:row_id, col_id_:col_id] <- mat
    col_id_ <- col_id + 1
  }
  for (i in 2:(n - 1)) {
    row_id_ <- 10 * (i - 1) + 1
    row_id <- 10 * (i - 1) + 10
    for (j in (i + 1):n) {
      col_id_ <- 10 * (j - 1) + 1
      col_id <- 10 * (j - 1) + 10
      cov[row_id_:row_id, col_id_:col_id] <- t(cov[1:10, row_id_:row_id]) %*%
        cov[1:10, col_id_:col_id]
    }
  }
  return(cov + t(cov))
}

test_mat <- create_cov_evals(10, 0.999)
any(eigen(test_mat)$values < 0)

image(test_mat)

t(test_mat[1:10, 11:20]) %*% test_mat[1:10, 21:30]


test_cov <- function(n_reps, n_views, alpha) {
  reps <- c()
  for (i in 1:n_reps) {
    # reps <- c(reps, min(eigen(create_cov_evals(8, ))$values))
    reps <- c(reps, any(eigen(create_cov_evals(n_views, alpha))$values < 0))
  }
  return(any(reps))
}

test_cov(3000, 8, 0.99) # 0.19
#'
test_cov(1000, 2, 0.99) # 1
test_cov(1000, 3, 0.50) # 0.501
test_cov(5000, 4, 0.3333) # 0.34
test_cov(5000, 5, 0.27) # 0.275
test_cov(5000, 6, 0.23) # 0.235
test_cov(5000, 7, 0.205) # 0.21
test_cov(5000, 8, 0.185) # 0.19
test_cov(5000, 9, 0.175) # 0.18
test_cov(5000, 10, 0.166) # 0.167
test_cov(500, 20, 0.111) # 0.12
test_cov(1000, 40, 0.07) # 0.08
 
n_views <- c(2:10, 20, 40)
alphas <- c(0.99, 0.50, 0.3333, 0.27, 0.23, 0.205, 0.185,
            0.175, 0.166, 0.111, 0.07)
plot(n_views, exp(-n_views), type = "l", xlab = "Number of Views",
  ylab = "Alpha", main = "Alpha vs Number of Views"
)
plot(n_views, alphas, type = "l", xlab = "Number of Views",
  ylab = "Alpha", main = "Alpha vs Number of Views"
)
plot(1/(n_views-1), alphas, type = "l", xlab = "Number of Views",
  ylab = "Alpha", main = "Alpha vs Number of Views"
)

(1/(n_views-1) - alphas )[11]
 plot(n_views[1:9], (1/(n_views-1) + 0.0072*n_views - alphas )[1:9])

# lines(n_views, alphas, col = "red")
legend("topright", legend = c("exp(-n_views)", "alphas"), col = c("black", "red"), lty = 1:2, cex = 0.8)

  type = "l", xlab = "Number of Views",
  ylab = "Alpha", main = "Alpha vs Number of Views"
)

exp(-(n_views)/(n_views - 1)) 

cov1 <- create_cov(10, 0.5)
any(eigen(create_cov(3, 0.5))$values < 0)

plot(n_views, (n_views-1)*exp(-((n_views)))/n_views, type = "l", xlab = "Number of Views",
  ylab = "Alpha", main = "Alpha vs Number of Views", ylim=c(0, 1)
)
lines(n_views, alphas, col = "red")

det(cov1)
det(create_cov_evals(19))

svd(create_cov_evals(2))$d
svd(cov1)$d


for (i in range(10)) {
  any(eigen(create_cov_evals(10, 0.5))$values < 0)
}







psd(cov1)
svd(cov1[1:10, 11:20])$d
diag(cov1[1:10, 11:20] %*% t(cov1[1:10, 11:20]))
cov1[1:10, 11:20] <- matrix(expm(A), 10, 10)
dim(expm(A))

reps <- c()
for (i in 1:5000) {
  reps <- c(reps, min(eigen(create_cov_evals(8, ))$values))
  # reps <- c(reps, any(eigen(create_cov_evals(10, 0.18))$values < 0))
}

any(reps < 0)
any(eigen(create_cov_evals(3, 0.6))$values < 0)
