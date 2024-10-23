#' ClusterCirc: Raw algorithm
#'
#' @description Is a help function that performs the ClusterCirc search algorithm
#'    within ClusterCirc-Data (cc_data) and ClusterCirc-Simu (cc_simu).
#'    All of the arguments are taken from cc_data or cc_simu. cc_raw cannot be
#'    performed on its own.
#'
#' @param angles Item angles to perform ClusterCirc on.
#' @param comm  Item communalities.
#' @param p Number of clusters.
#' @param m Number of variables.
#' @param w Vector with weights for the variables. Weights need to be in the same
#'    order as the variables in the data. Default = item communalities.
#' @param e Cluster weight (0 <= e <= 1) defining the importance of within-cluster
#'    proximity versus equal cluster spacing. Default is 1/p weighing all clusters
#'    equally. e = 0: Maximum importance of between-cluster spacing, within-cluster
#'    proximity is ignored. e = 1: Maximum importance of within-cluster proximity,
#'    between-cluster spacing is ignored.
#' @param q Precision index for the algorithm. Precision is higher for larger
#'   values. Default = 10.
#'
#' @return Returns item clusters with optimal circumplexity and ClusterCirc
#'   coefficients: Overall ClusterCirc results, coefficients for clusters and for items.
#'
#' @export


cc_raw <- function(angles, comm, p, m, w, e, q) {

  # ---------------------
  # ---- PREPARATION ----
  # ---------------------

  theta <- angles
  w_mn <- mean(w)

  # Sort theta and keep item number in ival (for later re-assignment)

  rk_th <- rank(theta)
  ival <- matrix(0, nrow = m, ncol = 3)

  for (i1 in 1:m) {
    for (i2 in 1:m) {
      if (rk_th[i1] == i2) {
        ival[i2, 1] <- i1
        ival[i2, 2] <- theta[i1]
        ival[i2, 3] <- w[i1]
      }
    }
  }

  # --------------------------------
  # ---- CLUSTERCIRC ALGORITHM ----
  # --------------------------------

  spacingw <- 361

  for (d in 0:round(360 * q / p)) {
    ci_h <- rep(0, m)
    c_m <- rep(0, p)
    c_no <- rep(0, p)
    c_min <- rep(0, p)
    c_max <- rep(0, p)
    c_rng <- rep(0, p)
    c_ang <- rep(0, p)

    # Check if item falls within the range of a cluster. Adjust for max. 360°.

    for (c in 1:p) {
      thmin <- (c - 1) * 360 / p + d / q
      thmax <- c * 360 / p + d / q

      for (i in 1:m) {
        if (thmin <= 360 & thmax <= 360) {
          if (ival[i, 2] >= thmin & ival[i, 2] < thmax) {
            ci_h[i] <- c
            c_m[c] <- c_m[c] + 1
          }
        }

        if (thmin <= 360 & thmax > 360) {
          thmax_n <- thmax - 360
          if ((ival[i, 2] >= thmin &
               ival[i, 2] <= 360) |
              (ival[i, 2] >= 0 & ival[i, 2] < thmax_n)) {
            ci_h[i] <- c
            c_m[c] <- c_m[c] + 1
          }
        }

        if (thmin > 360 & thmax > 360) {
          thmin_n <- thmin - 360
          thmax_n <- thmax - 360
          if (ival[i, 2] >= thmin_n & ival[i, 2] < thmax_n) {
            ci_h[i] <- c
            c_m[c] <- c_m[c] + 1
          }
        }
      }
    }

    # If c_m = 0: One cluster is empty. Prepare dismissal of solution.

    for (c in 1:p) {
      if (c_m[c] == 0) {
        c_m[c] <- 99
      }
    }

    # Compute cluster angle as the center between the outer items in cluster

    for (c in 1:p) {
      c_min[c] <- 361
      c_max[c] <- 0

      for (i in 1:m) {
        if (ci_h[i] == c & ival[i, 2] <= c_min[c]) {
          c_min[c] <- ival[i, 2]
        }
        if (ci_h[i] == c & ival[i, 2] >= c_max[c]) {
          c_max[c] <- ival[i, 2]
        }
      }

      if (c_max[c] - c_min[c] < 180) {
        c_ang[c] <- (c_max[c] + c_min[c]) / 2
        c_rng[c] <- c_max[c] - c_min[c]
      }

      # Special case: Cluster at approx. 0° could have a range of > 180.
      # Change item angles with help objects.

      if (c_max[c] - c_min[c] > 180) {
        ival_h <- ival[, 2]
        c_minh <- 361
        c_maxh <- 0

        for (i in 1:m) {
          if (ival[i, 2] > 180) {
            ival_h[i] <- ival[i, 2] - 360
          }
          if (ci_h[i] == c & ival_h[i] <= c_minh) {
            c_minh <- ival_h[i]
          }
          if (ci_h[i] == c & ival_h[i] >= c_maxh) {
            c_maxh <- ival_h[i]
          }
        }

        c_max[c] <- c_maxh
        c_min[c] <- 360 + c_minh
        c_rng[c] <- c_maxh - c_minh
        c_ang[c] <- (c_minh + c_maxh) / 2
      }
    }

    # Clusters and items need to be sorted before computing spacing indices

    c_rnk <- rank(c_ang)
    cval <-  matrix(0, nrow = p, ncol = 4)

    for (c1 in 1:p) {
      for (c2 in 1:p) {
        if (c_rnk[c2] == c1) {
          cval[c1, 1] <- c2
          cval[c1, 2] <- c_m[c2]
          cval[c1, 3] <- c_ang[c2]
          cval[c1, 4] <- c_rng[c2]
        }
      }
    }

    c_i <- rep(1, m)
    cvalh <- cval

    for (i in 1:m) {
      for (c in 1:p) {
        if (ci_h[i] == cval[c, 1]) {
          c_i[i] <- c
          cvalh[c, 1] <- c
        }
      }
    }

    # Help vectors for computation of distances:
    # Allow for negative angles (n) in first cluster

    ival_n <- cbind(c_i, ival)
    for (i in 1:m) {
      if (c_i[i] == 1 & ival[i, 2] > 180) {
        ival_n[i, 3] <- ival_n[i, 3] - 360
      }
    }

    rk_iang <- rank(ival_n[, 3])
    ival_h <- matrix(0, nrow = m, ncol = 4)
    for (i1 in 1:m) {
      for (i2 in 1:m) {
        if (rk_iang[i2] == i1) {
          ival_h[i1, ] = ival_n[i2, ]
        }
      }
    }

    # Compute spacing_w for each division

    ic_dis <- matrix(0, m, p)
    space <- 360 / p
    ic_dev <- matrix(0, m, p)
    ic_devp <- matrix(0, m, p)
    ic_dw <- matrix(0, m, p)
    ic_dwe <- matrix(0, m, p)

    for (i in 1:m) {
      for (c1 in 1:p) {
        c2 <- ival_h[i, 1]
        i_ang <- ival_h[i, 3]
        c_ang <- cvalh[c1, 3]
        i_w <- ival_h[i, 4]
        e_own <- e
        e_others <- (1-e)/(p-1)
        ic_dis[i, c1] <- i_ang - c_ang
        id_dis <- (c2 - 1) * space - (c1 - 1) * space
        ic_dev[i, c1] <- ic_dis[i, c1] - id_dis
        ic_devp[i, c1] <- ic_dev[i, c1] / space
        ic_dw[i, c1] <- ic_devp[i, c1] * sqrt(i_w)

        if (c1 == c2) {
          ic_dwe[i, c1] <- ic_dw[i, c1]*sqrt(e_own)
        }

        if (c1 != c2) {
          ic_dwe[i,c1] <- ic_dw[i, c1]*sqrt(e_others)
        }

      }
    }

    # Item spacing

    ispc_sq <- (apply(ic_devp ^ 2, 1, sum)) / p
    ispc <- sqrt(ispc_sq)

    # With weights for items (default = h_sq) and clusters (equal spacing assumption)
    # for spacing index (not interpretable on item level)

    w_mn <- mean(w)

    ispc_wsq <- (apply(ic_dwe ^ 2, 1, sum)) / w_mn
    ispc_w <- sqrt(ispc_wsq)

    # Overall spacing

    spc_sq <- (sum(ispc_sq)) / m
    spc <- sqrt(spc_sq)

    spc_wsq <- (sum(ispc_wsq)) / m
    spc_w <- sqrt(spc_wsq)

    # Dismiss partitions with empty clusters by making spc_w larger than
    # the initial spacingw (361). If this happens for all possible divisions,
    # the number of clusters is too large.

    for (c in 1:p) {
      if (c_m[c] == 99) {
        spc_w <- 50000
      }
    }

    if (spc_w < spacingw) {
      spacingw <- spc_w
      spacing <- spc
      items <- cbind(ival_h, ispc)
      clusters <- cvalh
      ic_dist <- ic_dis
    }

  }

  if (spacingw == 361) {
    print (
      "ClusterCirc could not finish, at least one of the clusters is empty.
       Try a smaller number of clusters or include more variables."
    )
  }

  if (spacingw < 361) {

    # Between-cluster spacing

    c_dis <- matrix(0, p, p)
    space <- 360 / p
    c_dev <- matrix(0, p, p)
    c_devp <- matrix(0, p, p)

    for (c1 in 1:p) {
      for (c2 in 1:p) {
        c_dis[c1, c2] <- clusters[c1, 3] - clusters[c2, 3]
        id_dis <- (c1 - 1) * space - (c2 - 1) * space
        c_dev[c1, c2] <- (c_dis[c1, c2] - id_dis)
        c_devp[c1, c2] <- c_dev[c1, c2] / space
      }
    }

    cbcs_psq <- (apply(c_devp ^ 2, 1, sum)) / (p - 1)
    cbcs_p <- sqrt(cbcs_psq)

    clusters <- cbind(clusters, cbcs_p)

    # Overall between-cluster spacing

    bcs_p <- sqrt(sum(cbcs_psq) / p)

    # Within-cluster proximity

    ic_disw <- rep(0, m)
    ic_disp <- rep(0, m)
    cwcp_sqp <- rep(0, p)
    cwcp_p <- rep(0, p)

    for (c in 1:p) {
      c_m <- clusters[c, 2]
      c_icdisp <- rep(0, m)
      for (i in 1:m) {
        if (items[i, 1] == c) {
          ic_disw[i] <- abs(ic_dist[i, c])
          ic_disp[i] <- ic_disw[i] / (360 / (2 * p))
          c_icdisp[i] <- ic_disp[i]
        }
      }

      cwcp_sqp[c] <- sum(c_icdisp ^ 2) / c_m
      cwcp_p[c] <- sqrt(cwcp_sqp[c])

    }

    # Overall within cluster spacing

    wcp_p <- sqrt(sum(cwcp_sqp) / p)

    # Final values

    clusters <- cbind(clusters, cwcp_p)
    items <- cbind(items, ic_disp)

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

    overall <- cbind(spacingw, spacing, bcs_p, wcp_p)

    # Parameters for cc_simu

    c_wrange <- clusters[, 4]
    wrange_c <- sum(c_wrange) / p
    for_cc_simu <- cbind(p, m, q, w_mn, wrange_c, e)

    output <- list(overall, clusters, items, for_cc_simu)

  }

}
