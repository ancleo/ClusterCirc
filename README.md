# ClusterCirc
This github repository features the R package "ClusterCirc", which finds
item clusters with optimal circumplex spacing in your data.

## Installing ClusterCirc in R:
Type the following commands in the R console to install and load ClusterCirc:

**install.packages("devtools")**  
**library(devtools)**    
**install_github("ancleo/ClusterCirc")**    
**library(ClusterCirc)**    

ClusterCirc depends on two external R packages: psych and knitr. If the
dependencies are not installed automatically, you can try the following:

**install_github("ancleo/ClusterCirc", dependencies = TRUE)**    
**library(ClusterCirc)**    

OR install the dependencies manually and re-install ClusterCirc by:

**install.packages("psych")**  
**library(psych)**  
**install.packages("knitr")**  
**library(knitr)**  
**install_github("ancleo/ClusterCirc")**  
**library(ClusterCirc)**  

To save installation time, the previous commands do not install the vignette, 
which contains a detailed description and demonstration of ClusterCirc functions.
If you want to install the vignette as well, type:  

**install_github("ancleo/ClusterCirc", build_vignettes = TRUE)**      
**library(ClusterCirc)**      

The description, vignette, and function documentation can be seen by  
**help(packages = "ClusterCirc")**    

## Downloading source code:

If you are interested in the source code of the ClusterCirc functions,
please download the current release of this repository. It contains all
files of the R package (Releases in the right-hand corner of this github website). 
Source code of the ClusterCirc functions can be found in the subfolder "R".

## Using and citing ClusterCirc:

The manuscript that presents ClusterCirc has been submitted to a peer-
reviewed journal. When using ClusterCirc, please cite the preprint version 
at the current stage of the publication process:

https://psyarxiv.com/yf37w/

(Note to reviewers: Please don't open the preprint version of the manuscript
to ensure double-blind review).

## Description of ClusterCirc

ClusterCirc is a clustering method designed for data with circular
structure. It can be used to find item clusters with optimal circumplex
spacing as an alternative to other clustering techniques like
conventional cluster analysis.

ClusterCirc is based on item angles. Input data can be item angles directly,
raw data, or item  loadings on two orthogonal factors or components from 
principal component analysis, exploratory or confirmatory factor analysis. 
If the data file contains raw data, item angles can be obtained by performing 
PCA (default) and converting the resulting loadings on two unrotated components 
into angles or by performing Browne's procedure CIRCUM on the data. 

ClusterCirc then sorts items into clusters that yield optimal circumplex
spacing. Optimal circumplex spacing for item clusters is given if clusters are
evenly distributed across the circle (equal spacing between clusters)
and if items are clustered closely around their cluster centroid
(minimal within-cluster spacing/item heterogeneity). Spacing coefficients 
are computed to assess circumplex spacing of items, clusters, and the overall data.
Range of all ClusterCirc coefficients: 0-1 (0 = perfect circumplex spacing).

### There are three functions for users:

1.  **cc_data:**  
    Main function. Sorts items of your dataset into
    clusters with optimal circumplex spacing. Spacing coefficients are
    computed for the suggested clustering. Depends on function cc_raw,
    which is included in the ClusterCirc package and automatically
    performed when cc_data is called.

    Usage on exemplary data with 3 clusters (p), 18 variables (n_var), items
    weighted by communalities (w_com = "TRUE", w), default precision index (q = 10):        
    **cc_data(file = data_ex, n_sample = 300, input = "PCA", p = 3,
        n_var = 18, w_com = TRUE, w, comm, e_def = TRUE, e, q = 10)**

3.  **cc_simu:**  
    Can be used to assess circumplex fit of the dataset.
    The function uses the specifications of the data and creates samples
    from a population with perfect circumplex spacing of clusters (default
    number of samples = 500). Results for the dataset (spacing coefficients 
    from cc_data) are compared to results from cc_simu to evaluate
    circumplexity in the data. cc_simu can only be used after performing
    cc_data.
	
	Usage for exemplary data (300 subjects in data):       
    **cc_simu(samples = 500, alpha = 1, input = "PCA")** 
	
4.  **cc_fix:** Computes ClusterCirc coefficients for user-defined item clusters
    without performing the ClusterCirc search algorithm to find item clusters.
    ClusterCirc-fix coefficients for user-defined item clusters can be compared
    to ClusterCirc coefficients for item clusters found by ClusterCirc-Data.
	
	Usage on exemplary data with Cluster1 = items 1 to 6, Cluster2 = items 7 to 10,
	Cluster3 = items 11 to 18 (limits):
    **cc_fix(file = data_ex, n_sample= 300, input = "PCA", limits = c(6,10,18),
       p = 3, n_var = 18, w_com = TRUE, comm, w, e_def = TRUE, e)** 

See function documentation in R for more detailed description and usage of functions:  
**?cc_data**    
**?cc_simu**    
**?cc_fix** 
