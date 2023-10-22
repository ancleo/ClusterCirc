#' Cluster-Circ: Raw algorithm
#'
#' @description Is a help function that performs the Cluster-Circ algorithm
#'    within Cluster-Circ Data (cc_data) and Cluster-Circ Simu (cc_simu).
#'    All of the arguments are taken from cc_data. cc_raw cannot be performed
#'    on its own.
#'
#' @param A Loading matrix to perform Cluster-Circ on.
#' @param p Number of clusters.
#' @param m Number of variables.
#' @param q Precision index for the algorithm. Precision is higher for larger
#'   values. Default = 10. Must be an integer > 0.
#'
#' @return Returns item clusters with optimal circumplexity and Cluster-Circ
#'   indices: Overall Cluster-Circ results, indices for clusters and for items.
#' @export
#'
#' @examples cc_raw(A, p = 3, m = 18, q = 10)

cc_raw <- function(A, p, m, q) {

  # ---------------------
  # ---- PREPARATION ----
  # ---------------------

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

  # Sort theta and keep item number in ival (for later re-assignment)

  rk_th <- rank(theta)
  ival <- matrix(0, nrow = m, ncol = 3)

  for (i1 in 1:m) {
    for (i2 in 1:m) {
      if (rk_th[i1] == i2) {
        ival[i2, 1] <- i1
        ival[i2, 2] <- theta[i1]
        ival[i2, 3] <- h_sq[i1]
      }
    }
  }

  # --------------------------------
  # ---- CLUSTER-CIRC ALGORITHM ----
  # --------------------------------

  spacingh <- 361

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

    # Compute spacing_h for each division

    ic_dis <- matrix(0, m, p)
    space <- 360 / p
    ic_dev <- matrix(0, m, p)
    ic_devp <- matrix(0, m, p)
    ic_dcom <- matrix(0, m, p)

    for (i in 1:m) {
      for (c1 in 1:p) {
        c2 <- ival_h[i, 1]
        i_ang <- ival_h[i, 3]
        c_ang <- cvalh[c1, 3]
        i_com <- ival_h[i, 4]
        ic_dis[i, c1] <- i_ang - c_ang
        id_dis <- (c2 - 1) * space - (c1 - 1) * space
        ic_dev[i, c1] <- ic_dis[i, c1] - id_dis
        ic_devp[i, c1] <- ic_dev[i, c1] / space
        ic_dcom[i, c1] <- ic_devp[i, c1] * sqrt(i_com)
      }
    }

    # Item spacing

    ispc_sq <- (apply(ic_devp ^ 2, 1, sum)) / p
    ispc <- sqrt(ispc_sq)

    # With communalities (h) for spacing index (not interpretable on item level)

    ispc_hsq <- (apply(ic_dcom ^ 2, 1, sum)) / (p * hsq_mn)
    ispc_h <- sqrt(ispc_hsq)

    # Overall spacing

    spc_sq <- (sum(ispc_sq)) / m
    spc <- sqrt(spc_sq)

    spc_hsq <- (sum(ispc_hsq)) / m
    spc_h <- sqrt(spc_hsq)

    # Dismiss partitions with empty clusters by making spc_h larger than
    # the initial spacingh (361). If this happens for all possible divisions,
    # the number of clusters is too large.

    for (c in 1:p) {
      if (c_m[c] == 99) {
        spc_h <- 50000
      }
    }

    if (spc_h < spacingh) {
      spacingh <- spc_h
      spacing <- spc
      items <- cbind(ival_h, ispc)
      clusters <- cvalh
      ic_dist <- ic_dis
    }

  }

  if (spacingh == 361) {
    print (
      "Cluster-Circ could not finish, at least one of the clusters is empty.
       Try a smaller number of clusters or include more variables."
    )
  }

  if (spacingh < 361) {

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

    # Final overall Cluster-Circ indices

    overall <- cbind(spacingh, spacing, bcs_p, wcp_p)

    # Parameters for cc_simu

    c_wrange <- clusters[, 4]
    wrange_c <- sum(c_wrange) / p
    for_cc_simu <- cbind(p, m, q, hsq_mn, wrange_c)

    output <- list(overall, clusters, items, for_cc_simu)

  }

}
