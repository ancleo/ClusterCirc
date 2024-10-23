#' ClusterCirc-Loadings
#'
#' @description Procedure that extracts item angles by trigonometric transformation
#'    of item loadings on two orthogonal axes.
#'
#' @param A Loading matrix on two orthogonal axes (e.g. from PCA or EFA).
#' @param m Number of variables.
#'
#' @return Returns item angles to be inserted into ClusterCirc.
#'
#' @export

cc_loadings <- function(A, m) {

  # Kaiser-normalization of loadings for simpler computation of angles

  h_sq <- apply(A ^ 2, 1, sum)
  h_rt <- sqrt(h_sq)
  hsq_mn <- mean(h_sq)
  Mh_rt <- matrix(diag(h_rt), ncol = m)
  A_k <- solve(Mh_rt) %*% A

  # Compute and sort theta: item angles in degrees
  # r = radians, p = positive

  A_pos <- abs(A_k)
  th_rp <- asin(A_pos[1:m, 1])
  th_r <- rep(0, m)
  theta <- rep(0, m)

  # Computation of angles depends on the quadrant (loadings positive/negative)

  for (i in 1:nrow(A_k)) {
    if (A_k[i, 1] >= 0 & A_k[i, 2] >= 0) {
      th_r[i] <- th_rp[i]
    }
    if (A_k[i, 1] >= 0 & A_k[i, 2] < 0) {
      th_r[i] <- pi - th_rp[i]
    }
    if (A_k[i, 1] < 0 & A_k[i, 2] < 0) {
      th_r[i] <- pi + th_rp[i]
    }
    if (A_k[i, 1] < 0 & A_k[i, 2] >= 0) {
      th_r[i] <- 2 * pi - th_rp[i]
    }
    theta[i] <- th_r[i] * 180 / pi
  }

  output <- list(theta,h_sq)

}
