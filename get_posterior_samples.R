#' Retrieve random samples of posterior coefficient estimates
#'
#' @param mod_list A list of brms regression model objects of class \code{brmsfit}.
#' @param L Integer indicating number of samples to pull
#'
#' @return a list of dataframes
#'
#' @author Xu Zhang
#'
#' @details The function loops through a list of *M* model objects, randomly samples coefficient estimates from the posterior, and compiles those samples in a single dataframe. Output is a list of *M* dataframes each with *n* columns of sampled estimates. Currently only works with \code{brmsfit} objects.
#' 
#' Additional details: this function changed the class() function which was used to detect model type to inherits () function. The inherits() function detects model families and hyper families. 
#' 
get_posterior_samples <- function(mod_list, n = 200L){
  
  post_samples <- lapply(mod_list, function(mod){

    ## ------------------------------------------------
    ## Detecting brms model
    ## ------------------------------------------------
    if (inherits(mod, "brmsfit")) {
      
      # extract all fixed-effect posterior draws
      post <- brms::posterior_samples(mod, pars = "^b_")
      
      # identify model family
      fam <- mod$family$family
      
      ### ------------------------------------------------
      ### Case 1: Binomial (Bernoulli) brms model
      ### ------------------------------------------------
      if (fam == "bernoulli") {
        
        sample_df <- post %>%
          dplyr::slice_sample(n = n) %>%
          t()
        
        return(sample_df)
      }
      
      ### ------------------------------------------------
      ### Case 2: Multinomial (Categorical) brms model
      ### ------------------------------------------------
      if (fam == "categorical") {
        
        # extract category names from parameter names
        cats <- unique(gsub("b_([^_]+)_.*", "\\1", names(post)))
        
        sample_by_cat <- lapply(cats, function(cat){
          
          # select parameters for this category
          post_cat <- post[, grep(paste0("^b_", cat, "_"), names(post)), drop = FALSE]
          
          # sample n posterior draws and transpose
          post_cat %>%
            dplyr::slice_sample(n = n) %>%
            t()
        })
        
        names(sample_by_cat) <- cats
        return(sample_by_cat)
      }
      
      ## ------------------------------------------------
      ## Unsupported brms family
      ## ------------------------------------------------
      stop("brms model family not supported.")
    }
    
    ## ------------------------------------------------
    ## Unsupported model class
    ## ------------------------------------------------
    stop("Model class not supported.")
  })
  
  names(post_samples) <- names(mod_list)
  return(post_samples)
}

