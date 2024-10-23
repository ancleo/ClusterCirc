#' ClusterCirc-Fix
#'
#' @description Computes ClusterCirc coefficients for user-defined item
#'    clusters. Items need to be sorted in the data file before analysis,
#'    e.g. cluster 1 = items 1 to 6, cluster 2 = items 7 to 10, etc.
#'    ClusterCirc-Fix coefficients for user-defined item clusters can be compared
#'    to ClusterCirc coefficients for item clusters found by ClusterCirc-Data.
#'
#' @param file File to be used for the analysis.
#' @param n_sample Sample size.
#' @param input Type of input to be inserted into ClusterCirc. Options are "PCA"
#'    (default), "Browne", "loadings", or "angles". "PCA" or "Browne" can be used
#'    if the file contains raw scores of the variables to obtain item angles from
#'    trigonometric transformation of two unrotated components (PCA) or from
#'    CIRCUM analysis (Browne). Insert "loadings" if the file contains loadings
#'    on two orthogonal axes from PCA or EFA or "angles" if it contains item angles.
#' @param limits Vector with the positions of the last item of each cluster.
#' @param p Number of clusters (minimum = 2).
#' @param n_var Number of variables.
#' @param w_com Is TRUE if communalities of the variables are used as weights.
#'    Is FALSE if user-defined weights should be used. Default = TRUE.
#' @param w Vector with weights for the variables. Weights need to be in the same
#'    order as the variables in the data and should be the same weights as used in
#'    cc_data for fair comparison of spc_w. Default = item communalities.
#' @param e_def Is TRUE if default values for cluster weights (importance of
#'    within-cluster proximity versus between-cluster spacing) are used.
#' @param e Cluster weight (0 <= e <= 1) defining the importance of within-cluster
#'    proximity versus equal cluster spacing. Default is 1/p weighing all clusters
#'    equally. e = 0: Maximum importance of between-cluster spacing, within-cluster
#'    proximity is ignored. e = 1: Maximum importance of within-cluster proximity,
#'    between-cluster spacing is ignored.
#' @param comm Vector with item communalities. Needs to be specified only if "input"
#'    is "angles". Otherwise, the argument is ignored, and item communalities are
#'    computed by the procedure.
#' @param w Vector with weights for the variables. Weights need to be in the same
#'    order as the variables in the data and should be the same weights as used in
#'    cc_data for fair comparison of spc_w. Default = item communalities.
#'
#' @return cc_fix returns overall ClusterCirc coefficients for user-defined item
#'    clusters: Overall ClusterCirc coefficients, coefficients for clusters and
#'    for items. Results can be compared to results for item clustering as
#'    suggested by cc_data.
#'
#' @export
#'
#' @examples cc_fix(file = data_ex, n_sample= 300, input = "PCA", limits = c(6,10,18),
#'    p = 3, n_var = 18, w_com = TRUE, comm, w, e_def = TRUE, e)

cc_fix <- function(file, n_sample, input = "PCA", limits, p, n_var, w_com = TRUE,
                   w, e_def = TRUE, e, comm) {

  # ---------------------
  # ---- PREPARATION ----
  # ---------------------

  if (input == "angles") {
    angles <- file
    comm <- comm
    cat("\n ClusterCirc-Fix is based on user-inserted item angles.")
  }

  if (input == "loadings") {
    cc_input <- cc_loadings(A = file, m = n_var)
    angles <- cc_input[[1]]
    comm <- cc_input[[2]]
    cat("\n Item angles are estimated by trigonometric transformation of")
    cat("\n user-inserted loadings on two components (circumplex axes).")
  }

  if (input == "PCA") {
    fit <- psych::principal(file, nfactors = 2, rotate = "none")
    fit[["loadings"]]
    A <- loadings(fit)
    class(A) <- "matrix"

    cc_input <- cc_loadings(A, m = n_var)
    angles <- cc_input[[1]]
    comm <- cc_input[[2]]

    cat("\n Item angles are estimated by trigonometric transformation of PCA")
    cat("\n loadings on two unrotated orthogonal components (circumplex axes).")
  }

  if (input == "Browne") {

    R_c <-as.matrix(round(cor(file),2))
    colnames(R_c) <- NULL
    rownames(R_c) <- NULL
    R_c[lower.tri(R_c)] <-0
    R_rd <- round(R_c,2)
    lmm = 10

    cc_browne <- tryCatch({
      circe_modfast(R_rd, v.names = seq(to = n_var), m = 1, N = n_sample, r = 1,
                    equal.com = FALSE, equal.ang = FALSE, mcsc = "unconstrained",
                    start.values = "PFA", file
      )
    }, error = function(msg) {
      return(matrix(NA, nrow = n_var * 3 + 1, ncol = 2))
    })

    if (any(is.na(cc_browne)) == TRUE) {
      cc_browne <- tryCatch({
        circe_det001(R_rd, v.names = seq(to = n_var), m = 1, N = n_sample, r = 1,
                     equal.com = FALSE, equal.ang = FALSE, mcsc = "unconstrained",
                     start.values = "PFA", file
        )
      }, error = function(msg) {
        return(matrix(NA, nrow = n_var * 3 + 1, ncol = 2))
      })
    }

    angles <- cc_browne[1:n_var,1]
    comm <- cc_browne[(n_var*2+2):(n_var*3+1),1]

    cat("\n Item angles are estimated by Browne's procedure: CIRCUM")
  }

  m <- n_var
  theta <- angles

  if (e_def == TRUE) {
    e <- 1/p
  }

  if (w_com == TRUE) {
    w <- comm
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

  c_no <- rep(1,m)
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
  ival <- matrix(0, nrow =m, ncol = 4)
  for (i1 in 1:m) {
    for (i2 in 1:m) {
      if (rk_iang[i2] == i1) {
        ival[i1, ] = ival_n[i2, ]
      }
    }
  }

  # -------------------------
  # -- CLUSTERCIRC INDICES --
  # -------------------------

  ic_dis <- matrix(0, m, p)
  space <- 360 / p
  ic_dev <- matrix(0, m, p)
  ic_devp <- matrix(0, m, p)
  ic_dw <- matrix(0, m, p)
  ic_dwe <- matrix(0, m, p)

  for (i in 1:m) {
    for (c1 in 1:p) {
      c2 <- ival[i, 1]
      i_ang <- ival[i, 3]
      c_ang <- cval[c1, 3]
      i_w <- ival[i, 4]
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

  cat("\n")
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


