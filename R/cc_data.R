#' ClusterCirc-Data
#'
#' @description Finds item clusters with optimal circumplex spacing in your data.
#'
#' @param file File to be used for the analysis.
#' @param n_sample Sample size.
#' @param input Type of input to be inserted into ClusterCirc. Options are "PCA"
#'    (default), "Browne", "loadings", or "angles". "PCA" or "Browne" can be used
#'    if the file contains raw scores of the variables to obtain item angles from
#'    trigonometric transformation of two unrotated components (PCA) or from
#'    CIRCUM analysis (Browne). Insert "loadings" if the file contains loadings
#'    on two orthogonal axes from PCA or EFA or "angles" if it contains item angles.
#' @param p Number of clusters (minimum = 2).
#' @param n_var Number of variables.
#' @param w_com Is TRUE if communalities of the variables are used as weights.
#'    Is FALSE if user-defined weights should be used. Default = TRUE.
#' @param w Vector with weights for the variables. Weights need to be in the same
#'    order as the variables in the data. Default = item communalities.
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
#' @param q Precision index for the algorithm. Precision is higher for larger
#'   values. Default = 10. Must be an integer > 0.
#'
#' @return Returns item clusters with optimal circumplexity and ClusterCirc
#'   coefficients: Overall ClusterCirc results, coefficients for clusters and for items.
#'
#' @export
#'
#' @examples cc_data(file = data_ex, n_sample = 300, input = "PCA", p = 3,
#'   n_var = 18, w_com = TRUE, w, comm, e_def = TRUE, e, q = 10)

cc_data <- function(file, n_sample, input = "PCA", p, n_var, w_com = TRUE, w,
                    e_def = TRUE, e, comm, q = 10) {

  if (input == "angles") {
    angles <- file
    comm <- comm
    cat("\n ClusterCirc is based on user-inserted item angles.")
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

  if (e_def == TRUE) {
    e <- 1/p
  }

  if (w_com == TRUE) {
    w <- comm
  }

  cc_results <- cc_raw(angles, comm, p, m = n_var, w, e, q)

  overall <- cc_results[[1]]
  clusters <- cc_results[[2]]
  items <- cc_results[[3]]
  for_cc_simu <- cc_results[[4]]

  for_cc_simu <- cbind(for_cc_simu, n_sample)

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
    "Weight (default: communality)",
    "Item-cluster spacing",
    "Distance to cluster center"
  )

  cc_overall <<- overall
  cc_clusters <<- clusters
  cc_items <<- items
  for_cc_simu <<- for_cc_simu

  cat("\n")
  cat("\n ==============================")
  cat("\n RESULTS CLUSTERCIRC-DATA")
  cat("\n ==============================", "\n")

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
  cat("The decision for the final item clusters is based on 'spacing (weighted)'.",
      "\n")
  cat("Weights are communalities of the items if not otherwise specified",
      "\n")
  cat("to account for (un-)reliability of measurement.",
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
