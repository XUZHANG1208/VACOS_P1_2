#' Calculate the second line of evidence for the VACOS method
#'
#' @param mod_list A list of brms model objects.
#' @param path Path in which to save the output as an R data file (\code{.rds}). If \code{NULL}, defaults to the current working directory. Set \code{path = FALSE} if you do not wish to save to file.
#' @param weight A numeric value indicating the size of the "effects" used for approximating the maximal reasonable distance. Default is 1.
#' @param scale How should the distance matrix be scaled? See details
#' @param overwrite Should the function overwrite data to location in \code{path}? Default is \code{'no'}, which will run the analysis if no file exists. If file in \code{path} exists, user with be prompted to set new path or allow file to be overwritten. Set to \code{'yes'} to automatically overwrite existing file, and \code{'reload'} to automatically reload existing file.
#' @param verbose Should messages be printed? Default is \code{FALSE}
#'
#' @author Xu Zhang
#'
#' @details The function loops through a list of model objects, extracts the coefficient estimates, and compiles them in a single dataframe.
#'
#' For scaling, there are four options. The default, \code{"abs"} (absolute), scales by a constant term based on the maximum reasonable distance, and values are bounded between 0 and 1 (see Szmrecsanyi et al. 2019). \code{"minmax"} uses minmax normalization, defined as
#'
#' \deqn{ x' = \frac{x - min(x)}{max(x) - min(x)}}{x' = (x - min(x))/(max(x) - min(x))}
#'
#' Minmax scaling bound values between 0 and 1. \code{"mean"} uses mean normalization, defined as

#' \deqn{ x' = \frac{x - mean(x)}{max(x) - min(x)}}{x' = (x - mean(x))/(max(x) - min(x))}
#'
#' If \code{scale = "none"} no scaling is applied.
#'
#' @return A \code{list} of length 3.
#' \describe{
#' \item{\code{coef.table}}{A dataframe of \emph{P} predictors by \emph{M} models, containing the pointwise estimated  coefficients (for \code{glm} and \code{glmer} models) or the mean posterior \beta estimates (for \code{brmsfit} models) for each predictor in each model.}
#' \item{\code{distance.matrix}}{An \emph{M} by \emph{M} distance matrix of class \code{dist}, derived from \code{coef.table}. Values are (normalized) Euclidean distances.}
#' \item{\code{similarity.scores}}{A dataframe of similarity scores derive from \code{distance.matrix}. See Szmrecsanyi et al. (2019) for details.}
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' data_list <- split(particle_verbs_short, particle_verbs_short$Variety, drop = TRUE)
#'
#' fmla <- Response ~ DirObjWordLength + DirObjDefiniteness + DirObjGivenness + DirObjConcreteness + DirObjThematicity + DirectionalPP + PrimeType + Semantics + Surprisal.P + Surprisal.V + Register
#'
#'brm_func <- function(x) brm(fmla, data = x, family = bernoulli/categorical)
#'
#' brm_list <- lapply(data_list, brm_func)
#' names(brm_list) <- names(data_list)
#'
#' line2 <- vacos_line2(brm_list, path = FALSE)
#' }
vacos_line2 <- function(mod_list, path = NULL, weight = 1, scale = c("abs", "mean", "minmax", "none"), overwrite = c("no",  "yes", "reload"), verbose = FALSE){
  
  overwrite <- match.arg(overwrite)
  
  if (is.null(path)) {
    path <- paste0(getwd(), "/vacos_line2_output_", format(Sys.time(), "%Y-%b-%d_%H-%M"), ".rds")
  }
  
  if(path == FALSE){
    output_list <- vector("list")
    create_coef_table <- function(mod_list) {
      coef_list <- lapply(mod_list, function(mod) {
        # ensure model is brms, including binomial and multinomial
        if (!inherits(mod, "brmsfit"))
          stop("All models must be brmsfit objects.")
        # extract posterior summary for all fixed effects
        post <- brms::posterior_summary(mod, pars = "^b_")
        # posterior means
        means <- post[, "Estimate"]
        # keep names
        names(means) <- rownames(post)
        return(means)
      })
      # combine into a table
      coef_tab <- do.call(cbind, coef_list)
      coef_tab <- as.data.frame(coef_tab)
      names(coef_tab) <- names(mod_list)
      return(coef_tab)
    }
    
    # Distance calculation
    raw_tab <- create_coef_table(mod_list)
    ## remove intercept if present
    coef_mat <- raw_tab[!grepl("Intercept", rownames(raw_tab)), ]
    ## transpose so models are rows
    dist_mat <- dist(t(coef_mat), method = "euclidean")
    
    
    if (match.arg(scale) == "abs"){
      # get the maximum reasonable distance
      dmy <- data.frame(a = sample(c(weight,-weight), size = nrow(coef_mat), replace = T))
      dmy$b <- -dmy$a # exact opposite of a
      maxD <- max(dist(t(dmy), "euclidean"))
      out_dist <- dist_mat/maxD
    } else if (match.arg(scale) == "minmax"){
      out_dist <- minmax(dist_mat)
    } else if (match.arg(scale) == "mean"){
      out_dist <- (dist_mat - mean(dist_mat))/(max(dist_mat) - min(dist_mat))
    } else {
      out_dist <- dist_mat
    }
    
    # Now normalize all distances to the maximum reasonable distance
    weighted_dist <- as.matrix(out_dist)
    diag(weighted_dist) <- NA # remove diagonals before calculating means
    means <- colMeans(weighted_dist, na.rm = T)
    sim_tab <- data.frame(Similarity = 1 - means)
    rownames(sim_tab) <- names(mod_list)
    
    # save normalized distances to output
    output_list[[2]] <- out_dist
    output_list[[3]] <- as.data.frame(sim_tab)
    
    names(output_list) <- c("coef.table", "distance.matrix", "similarity.scores")
  } else if(overwrite == "reload" & file.exists(path)){
    # reload from existing file
    if(verbose) message(paste("Loading existing file", path, "\nSet `overwrite = 'yes' or choose new path to calculate new values."))
    output_list <- readRDS(path)
  } else {
    output_list <- vector("list")
    raw_tab <- create_coef_table(mod_list) # call function to create varimp rankings
    output_list[[1]] <- raw_tab
    
    dist_mat <- dist(t(raw_tab[-1,]), method = "euclidean") # leave out the intercept
    
    if (match.arg(scale) == "abs"){
      # get the maximum reasonable distance
      dmy <- data.frame(a = sample(c(weight,-weight), size = nrow(raw_tab[-1,]), replace = T))
      dmy$b <- -dmy$a # exact opposite of a
      maxD <- max(dist(t(dmy), "euclidean"))
      out_dist <- dist_mat/maxD
    } else if (match.arg(scale) == "minmax"){
      out_dist <- minmax(dist_mat)
    } else if (match.arg(scale) == "mean"){
      out_dist <- (dist_mat - mean(dist_mat))/(max(dist_mat) - min(dist_mat))
    } else {
      out_dist <- dist_mat
    }
    
    # Now normalize all distances to the maximum reasonable distance
    weighted_dist <- as.matrix(out_dist)
    diag(weighted_dist) <- NA # remove diagonals before calculating means
    means <- colMeans(weighted_dist, na.rm = T)
    sim_tab <- data.frame(Similarity = 1 - means)
    rownames(sim_tab) <- names(mod_list)
    
    # save normalized distances to output
    output_list[[2]] <- out_dist
    output_list[[3]] <- as.data.frame(sim_tab)
    
    names(output_list) <- c("coef.table", "distance.matrix", "similarity.scores")
  }
  
  if(is.character(path)){
    if(overwrite == "yes"){
      if(file.exists(path) & verbose == TRUE) message("Existing file", path, "will be overwritten. Set overwrite = 'reload' to reload existing file.")
      saveRDS(output_list, file = path)
    } else if(overwrite == "no" & file.exists(path)) {
      msg <- paste("File", path, "already exists. Overwrite (y/n)?: ")
      over <- readline(prompt = msg)
      if(over == "y") {
        saveRDS(output_list, file = path)
      } else {
        new_path <- readline(prompt = "Please enter new file path:")
        saveRDS(output_list, file = new_path)
      }
    } else {
      saveRDS(output_list, file = path)
    }}
  
  return (output_list)
}
