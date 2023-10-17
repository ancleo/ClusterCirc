#' Cluster-Circ Simu

#' @description Simulates data with perfect circumplex clusters in the population
#'    using parameters from the empirical data (sample size, number of clusters,
#'    number of items, empirical within-cluster range, mean item communality).
#'    Performs Cluster-Circ on the simulated data for comparison with results
#'    from Cluster-Circ Data. Number of simulated samples: 500.
#'
#' @param n = Number of subjects in the sample.
#'
#' @return Returns Cluster-Circ indices for the population data and results for
#'    500 simulated samples (mean, SD, min, max). Empirical 'spacing with h²'
#'    from your data will be compared to the mean value from the simulated
#'    samples with perfect circumplexity to test whether circumplex fit of the
#'    empirical data is acceptable (using the cumulative probability of the
#'    standard normal distribution for p < .01, one-tailed).
#'
#' @export
#'
#' @examples cc_simu(n = 100)

cc_simu <- function(n) {

  cc_data_done <- exists("for_cc_simu")

  if (cc_data_done == FALSE) {
    cat("You must perform cc_data before performing cc_simu")

  } else {

    p <- for_cc_simu[1]
    m <- for_cc_simu[2]
    q <- for_cc_simu[3]
    h_sq <- for_cc_simu[4]
    h <- sqrt(rep(h_sq, m))
    c_wrange <- rep(for_cc_simu[5], p)
    n_simu <- 2
    n_pop <- n_simu * n
    mf <- m + 2
    # mf = Number of latent variables (m +2 circumplex factors)

    # ------------------------------
    # --- CREATE POPULATION DATA ---
    # ------------------------------

    # Normally distributed z-values for m error factors and 2 circumplex factors.
    # Perform PCA to orthogonalize factors and use regression factor scores (U, F).

    set.seed(1)
    Z <- replicate(mf, expr = rnorm(n_pop))
    fit <-
      psych::principal(Z,
                       nfactors = mf,
                       rotate = "none",
                       scores = TRUE)
    fscores <- fit[["scores"]]
    U <- fscores[, 1:m]
    F <- fscores[, (m + 1):mf]

    # The number of items per cluster mc in the population is chosen to be roughly
    # equal in all clusters. If m/p is not a natural number, some clusters have
    # fewer and some have more items:
    # mp_s = small mp, mp_l = large mp,  dec = decimals (for large and small mp)
    # mc = actual number of items in cluster c
    # p1 = number of clusters with small mc, p2 = number of clusters with large mc

    mp <- m / p
    mp_s <- trunc(mp)
    mp_l <- mp_s + 1
    dec_l <- mp - mp_s

    if (dec_l == 0) {
      mc <- rep(mp, p)
    }
    else if (dec_l > 0) {
      dec_s <- 1 - dec_l
      p1 <- p * dec_s
      p2 <- p * dec_l
      mc_s <- rep(mp_s, p1)
      mc_l <- rep(mp_l, p2)
      mc <- cbind(mc_s, mc_l)
    }

    # Introduce circumplexity: Translate perfect circumplex angles into population
    # loadings. Use reg factor scores and make data with var = F*T(A) + U*T(D).

    cl_dis <- 360 / p
    ones <- rep(1, p)
    w_dis <- c_wrange %/% (mc - ones)

    for (c in 1:p) {
      th <- rep(mc[c], 0)
      for (i in 1:mc[c]) {
        th[i] <- cl_dis * (c - 1) + (i - (mc[c] + 1) / 2) * w_dis[c]
      }
      if (c == 1) {
        theta_d <- th
      }
      else if (c > 1) {
        theta_d <- c(theta_d, th)
      }
    }

    for (i in 1:m) {
      if (theta_d[i] < 0) {
        theta_d[i] <- 360 + theta_d[i]
      }
    }

    theta <- theta_d * (pi / 180)

    # Angles > 2*pi (> 360°) need to be converted

    for (i in 1:m) {
      if (theta[i] > 2 * pi) {
        theta[i] <- theta[i] - 2 * pi
      }
    }

    # Conversion with sin/cos works only for theta >= 0 and >= pi/2 (90°).
    # Help vector theta_h for compuation of population loadings.

    theta_h <- theta
    for (i in 1:m) {
      if (theta[i] > pi / 2 & theta[i] <= pi) {
        theta_h[i] <- pi - theta[i]
      }
      if (theta[i] > pi & theta[i] <= 3 / 2 * pi) {
        theta_h[i] <- theta[i] - pi
      }
      if (theta[i] > 3 / 2 * pi & theta[i] <= 2 * pi) {
        theta_h[i] <- 2 * pi - theta[i]
      }
    }

    A = matrix(0, nrow = m, ncol = 2)

    for (i in 1:m) {
      if (theta[i] >= 0 & theta[i] <= pi / 2) {
        A[i, 1] <- h[i] * sin(theta_h[i])
        A[i, 2] <- h[i] * cos(theta_h[i])
      }
      if (theta[i] > pi / 2 & theta[i] <= pi) {
        A[i, 1] <- h[i] * sin(theta_h[i])
        A[i, 2] <- -(h[i] * cos(theta_h[i]))
      }
      if (theta[i] > pi & theta[i] <= 3 / 2 * pi) {
        A[i, 1] <- -(h[i] * sin(theta_h[i]))
        A[i, 2] <- -(h[i] * cos(theta_h[i]))
      }
      if (theta[i] > 3 / 2 * pi & theta[i] <= 2 * pi) {
        A[i, 1] <- -(h[i] * sin(theta_h[i]))
        A[i, 2] <- h[i] * cos(theta_h[i])
      }
    }

    # Get error loadings D for var = F*T(A) + U*T(D)

    ones <- rep(1, m)
    ones <- matrix(diag(ones), ncol = m)
    AtA <- A %*% t(A)
    AtA <- diag(diag(AtA))
    D <- sqrt(ones - AtA)

    # Compute variables (simulated values for each subject in the population):

    var <- F %*% t(A) + U %*% t(D)

    # ------------------------------
    # --- PCA ON POPULATION DATA ---
    # ------------------------------

    fit <- psych::principal(var, nfactors = 2, rotate = "none")
    fit[["loadings"]]
    A_pop <- loadings(fit)
    class(A_pop) <- "matrix"
    A_pop

    # -------------------------------------------
    # --- CLUSTER-CIRC ON POPULATION LOADINGS ---
    # -------------------------------------------

    cc_pop <- cc_raw(A_pop, p, m, q)

    overall_pop <- cc_pop[[1]]
    clusters_pop <- cc_pop[[2]]
    items_pop <- cc_pop[[3]]

    # ------------------------
    # --- SORTING CORRECT? ---
    # ------------------------

    # clus_id = ideal cluster sorting

    for (c in 1:p) {
      c_id <- rep(c, mc[c])
      if (c == 1) {
        clus_id <- c_id
      }
      else if (c > 1) {
        clus_id <- c(clus_id, c_id)
      }
    }

    # Each entry in 'sort_it' indicates whether items are clustered together
    # correctly (1) or wrongly (0). A single zero entry indicates that the
    # intended structure has not been found.
    # Position in clus_id, items[i, 2] = item number.
    # items[i, 1] is the cluster to which the item has been sorted.

    sort_it <- matrix(0, nrow = m, ncol = m)

    for (i1 in 1:m) {
      for (i2 in 1:m) {
        if (i1 != i2) {
          if (clus_id[items_pop[i1, 2]] == clus_id[items_pop[i2, 2]] &
              items_pop[i1, 1] == items_pop[i2, 1]) {
            sort_it[i1, i2] <- 1
          }
          if (clus_id[items_pop[i1, 2]] == clus_id[items_pop[i2, 2]] &
              items_pop[i1, 1] != items_pop[i2, 1]) {
            sort_it[i1, i2] <- 0
          }
          if (clus_id[items_pop[i1, 2]] != clus_id[items_pop[i2, 2]] &
              items_pop[i1, 1] == items_pop[i2, 1]) {
            sort_it[i1, i2] <- 0
          }
          if (clus_id[items_pop[i1, 2]] != clus_id[items_pop[i2, 2]] &
              items_pop[i1, 1] != items_pop[i2, 1]) {
            sort_it[i1, i2] <- 1
          }
        }
      }
    }

    if (sum(rowSums(sort_it)) == m * (m - 1)) {
      sortcorr_pop <- 1
    }
    else if (sum(rowSums(sort_it)) < m * (m - 1)) {
      sortcorr_pop <- 0
    }

    as.data.frame(overall_pop)
    colnames(overall_pop) = c(
      "Spacing (with h²)",
      "Spacing",
      "Between-cluster spacing",
      "Within-cluster proximity"
    )

    as.data.frame(clusters_pop)
    colnames(clusters_pop) = c(
      "Cluster",
      "Items",
      "Angle",
      "Range",
      "Between-cluster spacing",
      "Within-cluster proximity"
    )

    as.data.frame(items_pop)
    colnames(items_pop) = c(
      "Cluster",
      "Item",
      "Angle",
      "Communality",
      "Item-cluster spacing",
      "Distance to cluster center"
    )

    # --------------------------------
    # --- PCA ON SIMULATED SAMPLES ---
    # --------------------------------

    A_samp <- matrix(0, nrow = n_simu * m, ncol = 2)

    for (i_simu in 1:n_simu) {
      var_s <- var[((i_simu - 1) * n + 1):(i_simu * n), ]
      fit <- psych::principal(var_s, nfactors = 2, rotate = "none")
      fit[["loadings"]]
      A_s <- loadings(fit)
      class(A_s) <- "matrix"
      A_s

      # -------------------------------
      # --- CLUSTER-CIRC ON SAMPLES ---
      # -------------------------------

      cc_s <- cc_raw(A_s, p, m, q)

      overall_s <- cc_s[[1]]
      clusters_s <- cc_s[[2]]
      items_s <- cc_s[[3]]

      # ------------------------
      # --- SORTING CORRECT? ---
      # ------------------------

      for (c in 1:p) {
        c_id <- rep(c, mc[c])
        if (c == 1) {
          clus_id <- c_id
        }
        else if (c > 1) {
          clus_id <- c(clus_id, c_id)
        }
      }

      sort_it <- matrix(0, nrow = m, ncol = m)

      for (i1 in 1:m) {
        for (i2 in 1:m) {
          if (i1 != i2) {
            if (clus_id[items_s[i1, 2]] == clus_id[items_s[i2, 2]] &
                items_s[i1, 1] == items_s[i2, 1]) {
              sort_it[i1, i2] <- 1
            }
            if (clus_id[items_s[i1, 2]] == clus_id[items_s[i2, 2]] &
                items_s[i1, 1] != items_s[i2, 1]) {
              sort_it[i1, i2] <- 0
            }
            if (clus_id[items_s[i1, 2]] != clus_id[items_s[i2, 2]] &
                items_s[i1, 1] == items_s[i2, 1]) {
              sort_it[i1, i2] <- 0
            }
            if (clus_id[items_s[i1, 2]] != clus_id[items_s[i2, 2]] &
                items_s[i1, 1] != items_s[i2, 1]) {
              sort_it[i1, i2] <- 1
            }
          }
        }
      }

      if (sum(rowSums(sort_it)) == m * (m - 1)) {
        sortcorr <- 1
      }
      else if (sum(rowSums(sort_it)) < m * (m - 1)) {
        sortcorr <- 0
      }

      # -------------------------------
      # --- COMBINE FOR ALL SAMPLES ---
      # -------------------------------

      if (i_simu == 1) {
        ovrl_all <- overall_s
        sort_all <- sortcorr
      }
      if (i_simu > 1) {
        ovrl_all <- rbind(ovrl_all, overall_s)
        sort_all <- rbind(sort_all, sortcorr)
      }

    }

    mean_simu <- colMeans(ovrl_all)
    sd_simu <- apply(ovrl_all, 2, sd)
    min_simu <- apply(ovrl_all, 2, min)
    max_simu <- apply(ovrl_all, 2, max)

    overall_simu <- rbind(mean_simu, sd_simu, min_simu, max_simu)

    accuracy <- sum(sort_all) / n_simu * 100

    simu_par <- matrix(c(n_simu, n, m, p, q), nrow = 1, ncol = 5)

    as.data.frame(overall_simu)
    colnames(overall_simu) = c(
      "Spacing (with h²)",
      "Spacing",
      "Between-cluster spacing",
      "Within-cluster proximity"
    )
    rownames(overall_simu) = c("Mean", "SD", "Minimum", "Maximum")

    as.data.frame(simu_par)
    colnames(simu_par) = c("Number of samples",
                           "Sample size",
                           "Items",
                           "Clusters",
                           "Precision index")

    spch_dat <- cc_overall[1]
    spch_sim <- mean_simu[1]
    spch_sd = sd_simu[1]
    cutoff <- spch_sim + 2.33 * spch_sd
    comp <- cutoff - spch_dat

    cat("\n")
    cat("============================", "\n")
    cat("RESULTS CLUSTER-CIRC SIMU", "\n")
    cat("============================", "\n")

    cat("\n")
    cat("PARAMETERS OF THE SIMULATION")
    print(knitr::kable(simu_par, "simple", digits = 3))

    cat("\n")
    cat("============================", "\n")
    cat("POPULATION (SIMULATED)", "\n")
    cat("============================", "\n")

    cat("\n")
    cat("The following tables show results for a simulated population with", "\n")
    cat("perfect circumplex spacing of clusters, adapted for the specifications", "\n")
    cat("of the empirical data: Sample size, number of clusters, number of items,", "\n")
    cat("empirical within-cluster range, mean item communality.", "\n")

    cat("\n")
    cat("OVERALL - POPULATION")
    print(knitr::kable(overall_pop, "simple", digits = 3))

    cat("\n")
    cat("CLUSTERS - POPULATION")
    print(knitr::kable(clusters_pop, "simple", digits = 3))

    cat("\n")
    cat("ITEMS - POPULATION")
    print(knitr::kable(items_pop, "simple", digits = 3))

    if (sortcorr_pop == 1) {
      cat("\n")
      cat("Cluster-Circ found the intended circumplex clusters in the population.",
          "\n")
    } else if (sortcorr_pop == 0) {
      cat("\n")
      cat("Cluster-Circ did not find the intended circumplex clusters in the
          population.","\n")
    }

    cat("\n")
    cat("==================================", "\n")
    cat("SAMPLES OF THE SIMULATION STUDY", "\n")
    cat("==================================", "\n")

    cat("\n")
    cat("The following tables show results for 500 simulated samples from the population",
        "\n")
    cat("with perfect circumplex clusters.", "\n")

    cat("\n")
    cat("ACCURACY:", accuracy, "%", "\n")
    cat("of the simulated samples sorted the items according to the intended circumplex clusters.",
        "\n")

    cat("\n")
    cat("OVERALL RESULTS OF THE SIMULATED SAMPLES")
    print(knitr::kable(overall_simu, "simple", digits = 3))

    cat("\n")
    cat("Range of all Cluster-Circ coefficients: 0-1 (0 = perfect circumplex spacing).",
        "\n")

    cat("\n")
    cat("Spacing with h² in empirical data:", spch_dat, "\n")

    cat("\n")
    cat("Recommendation: Circumplex fit of the empirical data is acceptable if 'spacing",
        "\n")
    cat("with h²' in the empirical data is not larger than mean 'spacing with h²' + 2.33 SD",
        "\n")
    cat("from the simulated samples in Cluster-Circ Simu (corresponding to the cumulative",
        "\n")
    cat("probability of the standard normal distribution for p < .01, one-tailed).",
        "\n")

    if (comp >= 0) {
      cat("\n")
      cat("Here: Empirical 'spacing with h²' is within mean-spc_h + 2.33 SD from",
          "\n")
      cat("Cluster-Circ Simu -> Circumplex fit acceptable.", "\n")
    }
    else if (comp < 0) {
      cat("\n")
      cat("Here: Empirical 'spacing with h²' is larger than mean-spc_h + 2.33 SD",
          "\n")
      cat("from Cluster-Circ Simu -> Low circumplex fit.", "\n")
    }

    cat("\n")
    cat("The manuscript that presents Cluster-Circ has been submitted to a peer-",
        "\n")
    cat ("reviewed journal. When using Cluster-Circ, please cite the preprint version",
         "\n")
    cat ("at the current stage of the publication process:", "\n")
    cat("\n")
    cat("Weide et al. (2022): Cluster-Circ: Finding item clusters for circumplex",
        "\n")
    cat ("instruments. PsyArxiv (preprint). https://psyarxiv.com/yf37w/")

  }

}
