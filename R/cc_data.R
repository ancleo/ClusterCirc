#' ClusterCirc-Data
#'
#' @description Finds item clusters with optimal circumplex spacing in your data.
#'    If the file contains raw scores (type = "scores"), PCA without rotation
#'    will be performed on the data before running ClusterCirc.
#'
#' @param file File to be used for the analysis.
#' @param type Fill in "scores" if the file contains raw scores of variables,
#'    "loadings" if the file contains loadings from PCA or EFA.
#'    Default is "scores".
#' @param p Number of clusters (minimum = 2).
#' @param m Number of variables.
#' @param w_com Is "TRUE" if communalities of the variables are used as weights.
#'    Is "FALSE" if user-defined weights should be used. Default = "TRUE".
#' @param w Vector with weights for the variables. Weights need to be in the same
#'    order as the variables in the data. Default = item communalities.
#' @param q Precision index for the algorithm. Precision is higher for larger
#'   values. Default = 10. Must be an integer > 0.
#'
#' @return Returns item clusters with optimal circumplexity and ClusterCirc
#'   coefficients: Overall ClusterCirc results, coefficients for clusters and for items.
#' @export
#'
#' @examples cc_data(file = data_ex, type = "scores", p = 3, m = 18, w_com = "TRUE", w, q = 10)
#'
cc_data <- function(file, type = "scores", p, m, w_com = "TRUE", w, q = 10) {

  if (type == "scores") {
    fit <- psych::principal(file, nfactors = 2, rotate = "none")
    fit[["loadings"]]
    A <- loadings(fit)
    class(A) <- "matrix"
  }

  if (type == "loadings") {
    A <- file
  }

  cc_results <- cc_raw(A, p, m, w_com, w, q)

  overall <- cc_results[[1]]
  clusters <- cc_results[[2]]
  items <- cc_results[[3]]
  for_cc_simu <- cc_results[[4]]

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

  cat("\n ==============================")
  cat("\n RESULTS CLUSTERCIRC-DATA ")
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
