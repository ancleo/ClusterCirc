#' ClusterCirc-fix
#'
#' @description Computes ClusterCirc coefficients for user-defined item
#'    clusters. Items need to be sorted in the data file before analysis,
#'    e.g. cluster 1 = items 1 to 6, cluster 2 = items 7 to 10, etc.
#'    ClusterCirc-fix coefficients for user-defined item clusters can be compared
#'    to ClusterCirc coefficients for item clusters found by ClusterCirc-Data.
#'
#' @param file File to be used for the analysis.
#' @param type Is "scores" if the file contains raw scores of variables,
#'    "loadings" if the file contains loadings from PCA or EFA.
#'    Default is "scores".
#' @param limits Vector with the positions of the last item of each cluster,
#'    e.g. c(6,10,18) for 3 clusters with C1 = i1-i6, C2 = i7-i10, C3 = i11-i18.
#' @param p Number of clusters.
#' @param m Number of variables.
#' @param w_com Is "TRUE" if communalities of the variables are used as weights.
#'    Is "FALSE" if user-defined weights should be used. Default = "TRUE".
#' @param w Vector with weights for the variables. Weights need to be in the same
#'    order as the variables in the data and should be the same weights as used in
#'    cc_data for fair comparison of spc_w. Default = item communalities.
#'
#' @return cc_fix returns overall ClusterCirc coefficients for user-defined item
#'    clusters: Overall ClusterCirc coefficients, coefficients for clusters and
#'    for items. Results can be compared to results for item clustering as
#'    suggested by cc_data.
#' @export
#'
#' @examples cc_fix(data_ex, type = "scores", limits = c(6,10,18), p = 3, m = 18, w_com = "TRUE", w)

cc_fix <- function(file, type = "scores", limits, p, m, w_com = "TRUE", w) {

  # ---------------------
  # ---- PREPARATION ----
  # ---------------------

  if (type == "scores") {
    fit <- psych::principal(file, nfactors = 2, rotate = "none")
    fit[["loadings"]]
    A <- loadings(fit)
    class(A) <- "matrix"
  }

  if (type == "loadings") {
    A <- file
  }

  # Kaiser-normalization of loadings for simpler computation of angles

  h_sq <- apply(A ^ 2, 1, sum)
  h_rt <- sqrt(h_sq)
  hsq_mn <- mean(h_sq)
  Mh_rt <- matrix(diag(h_rt), ncol = m)
  A_k <- solve(Mh_rt) %*% A

  if (w_com == "TRUE")
    w <- h_sq
  w_mn <- hsq_mn

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

  # Fixed item-cluster assignment according to 'limits'

  item_no <- rep(1:m)
  c_m <- rep(0, p)

  for (c in 1:p) {
    if (c == 1) {
      c_m[c] <- limits[c]
      ci <- rep(c, limits[c])
      ci_v <- ci
    }
    if (c > 1) {
      c_m[c] <- limits[c] - limits[c - 1]
      ci <- rep(c, limits[c] - limits[c - 1])
      ci_v <- c(ci_v, ci)
    }
  }

  # Compute cluster angle as the center between the outer items in cluster

  c_min <- rep(0, p)
  c_max <- rep(0, p)
  c_rng <- rep(0, p)
  c_ang <- rep(0, p)

  for (c in 1:p) {
    c_min[c] <- 361
    c_max[c] <- 0

    for (i in 1:m) {
      if (ci_v[i] == c & theta[i] <= c_min[c]) {
        c_min[c] <- theta[i]
      }
      if (ci_v[i] == c & theta[i] >= c_max[c]) {
        c_max[c] <- theta[i]
      }
    }

    if (c_max[c] - c_min[c] < 180) {
      c_ang[c] <- (c_max[c] + c_min[c]) / 2
      c_rng[c] <- c_max[c] - c_min[c]
    }

    # Special case: Cluster at approx. 0° could have a range of > 180.
    # Change item angles with help objects.

    if (c_max[c] - c_min[c] > 180) {
      theta_h <- theta
      c_minh <- 361
      c_maxh <- 0

      for (i in 1:m) {
        if (theta[i] > 180) {
          theta_h[i] <- theta[i] - 360
        }
        if (ci_v[i] == c & theta_h[i] <= c_minh) {
          c_minh <- theta_h[i]
        }
        if (ci_v[i] == c & theta_h[i] >= c_maxh) {
          c_maxh <- theta_h[i]
        }
      }

      c_max[c] <- c_maxh
      c_min[c] <- 360 + c_minh
      c_rng[c] <- c_maxh - c_minh
      c_ang[c] <- (c_minh + c_maxh) / 2
    }
  }

  # Clusters and items need to be sorted before computing spacing indices

  ival_h <- cbind(item_no, theta, w)
  cval_h <-  matrix(0, nrow = p, ncol = 4)

  c_rnk <- rank(c_ang)

  for (c1 in 1:p) {
    for (c2 in 1:p) {
      if (c_rnk[c2] == c1) {
        cval_h[c1, 1] <- c2
        cval_h[c1, 2] <- c_m[c2]
        cval_h[c1, 3] <- c_ang[c2]
        cval_h[c1, 4] <- c_rng[c2]
      }
    }
  }

  c_no <- rep(1, m)
  cval <- cval_h

  for (i in 1:m) {
    for (c in 1:p) {
      if (ci_v[i] == cval_h[c, 1]) {
        c_no[i] <- c
        cval[c, 1] <- c
      }
    }
  }

  ival_n <- cbind(c_no, ival_h)

  # Help vectors for computation of distances:
  # Allow for negative angles (n) in first cluster

  for (i in 1:m) {
    if (c_no[i] == 1 & ival_h[i, 2] > 180) {
      ival_n[i, 3] <- ival_n[i, 3] - 360
    }
  }

  rk_iang <- rank(ival_n[, 3])
  ival <- matrix(0, nrow = m, ncol = 4)
  for (i1 in 1:m) {
    for (i2 in 1:m) {
      if (rk_iang[i2] == i1) {
        ival[i1, ] = ival_n[i2, ]
      }
    }
  }

  # ------------------------------
  # ---- CLUSTERCIRC INDICES ----
  # ------------------------------

  ic_dis <- matrix(0, m, p)
  space <- 360 / p
  ic_dev <- matrix(0, m, p)
  ic_devp <- matrix(0, m, p)
  ic_dw <- matrix(0, m, p)

  for (i in 1:m) {
    for (c1 in 1:p) {
      c2 <- ival[i, 1]
      i_ang <- ival[i, 3]
      c_ang <- cval[c1, 3]
      i_w <- ival[i, 4]
      ic_dis[i, c1] <- i_ang - c_ang
      id_dis <- (c2 - 1) * space - (c1 - 1) * space
      ic_dev[i, c1] <- ic_dis[i, c1] - id_dis
      ic_devp[i, c1] <- ic_dev[i, c1] / space
      ic_dw[i, c1] <- ic_devp[i, c1] * sqrt(i_w)
    }
  }

  # Item spacing

  ispc_sq <- (apply(ic_devp ^ 2, 1, sum)) / p
  ispc <- sqrt(ispc_sq)

  # With communalities (h) for spacing index (not interpretable on item level)


  ispc_wsq <- (apply(ic_dw ^ 2, 1, sum)) / (p * w_mn)
  ispc_w <- sqrt(ispc_wsq)

  # Overall spacing

  spc_sq <- (sum(ispc_sq)) / m
  spc <- sqrt(spc_sq)

  spc_wsq <- (sum(ispc_wsq)) / m
  spc_w <- sqrt(spc_wsq)

  # Between-cluster spacing

  c_dis <- matrix(0, p, p)
  space <- 360 / p
  c_dev <- matrix(0, p, p)
  c_devp <- matrix(0, p, p)

  for (c1 in 1:p) {
    for (c2 in 1:p) {
      c_dis[c1, c2] <- cval[c1, 3] - cval[c2, 3]
      id_dis <- (c1 - 1) * space - (c2 - 1) * space
      c_dev[c1, c2] <- (c_dis[c1, c2] - id_dis)
      c_devp[c1, c2] <- c_dev[c1, c2] / space
    }
  }

  cbcs_psq <- (apply(c_devp ^ 2, 1, sum)) / (p - 1)
  cbcs_p <- sqrt(cbcs_psq)

  # Overall between-cluster spacing

  bcs_p <- sqrt(sum(cbcs_psq) / p)

  # Within-cluster proximity

  ic_disw <- rep(0, m)
  ic_disp <- rep(0, m)
  cwcp_sqp <- rep(0, p)
  cwcp_p <- rep(0, p)

  for (c in 1:p) {
    c_m <- cval[c, 2]
    c_icdisp <- rep(0, m)
    for (i in 1:m) {
      if (ival[i, 1] == c) {
        ic_disw[i] <- abs(ic_dis[i, c])
        ic_disp[i] <- ic_disw[i] / (360 / (2 * p))
        c_icdisp[i] <- ic_disp[i]
      }
    }

    cwcp_sqp[c] <- sum(c_icdisp ^ 2) / c_m
    cwcp_p[c] <- sqrt(cwcp_sqp[c])

  }

  # Overall within-cluster proximity:

  wcp_p <- sqrt(sum(cwcp_sqp) / p)

  # Final values

  items <- cbind(ival, ispc, ic_disp)
  clusters <- cbind(cval, cbcs_p, cwcp_p)

  # Recompute negative angles

  for (i in 1:m) {
    if (items[i, 3] < 0) {
      items[i, 3] <- 360 + items[i, 3]
    }
  }

  for (c in 1:p) {
    if (clusters[c, 3] < 0) {
      clusters[c, 3] <- 360 + clusters[c, 3]
    }
  }

  # Final overall ClusterCirc indices

  overall <- cbind(spc_w, spc, bcs_p, wcp_p)

  as.data.frame(overall)
  colnames(overall) = c(
    "Spacing (weighted)",
    "Spacing",
    "Between-cluster spacing",
    "Within-cluster proximity"
  )

  as.data.frame(clusters)
  colnames(clusters) = c(
    "Cluster",
    "Items",
    "Angle",
    "Range",
    "Between-cluster spacing",
    "Within-cluster proximity"
  )

  as.data.frame(items)
  colnames(items) = c(
    "Cluster",
    "Item",
    "Angle",
    "Weight (default = communality)",
    "Item-cluster spacing",
    "Distance to cluster center"
  )

  cc_fix_overall <<- overall
  cc_fix_clusters <<- clusters
  cc_fix_items <<- items

  cat("\n ============================")
  cat("\n RESULTS CLUSTERCIRC-FIX ")
  cat("\n ============================", "\n")

  cat("\n OVERALL")
  print(knitr::kable(overall, "simple", digits = 3))
  cat("\n CLUSTERS")
  print(knitr::kable(clusters, "simple", digits = 3))
  cat("\n ITEMS")
  print(knitr::kable(items, "simple", digits = 3))

  cat("\n")
  cat("Range of all ClusterCirc coefficients: 0-1 (0 = perfect circumplex spacing).",
      "\n")
  cat("\n")
  cat("Item weights are communalities of the items if not otherwise specified.",
      "\n")
  cat("\n")
  cat("The manuscript that presents ClusterCirc has been submitted to a peer-",
      "\n")
  cat ("reviewed journal. When using ClusterCirc, please cite the preprint version",
       "\n")
  cat ("at the current stage of the publication process:", "\n")
  cat("\n")
  cat("Note to reviewers: Please don't open the preprint article during blind peer-review.",
      "\n")
  cat ("https://psyarxiv.com/yf37w/")

}
