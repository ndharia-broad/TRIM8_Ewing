
## Publication theme ggplot
theme_Publication <- function(base_size=14, base_family="Helvetica") {
  theme_foundation(base_size=base_size, base_family=base_family) +
    theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(), 
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          # panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.justification = 'top',
          legend.title = element_text(face="italic"),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold")
    )
}

## Pubmed search
pubmed_search <- function(search_term1, search_term2) {
  search_term <- paste0(search_term1, " ", 
                        gsub("t cell", "T cell",
                             gsub("b cell", "B cell",
                                  gsub("HER2pos", "HER2+", 
                                       gsub("HER2neg", "HER2-", 
                                            gsub("ERneg", "ER-",
                                                 gsub("ERpos", "ER+", 
                                                      gsub("_", " ", search_term2))))))), 
                        "")
  search_info <- NULL
  sleep <- 0.5
  attempt <- 1
  
  while(is.null(search_info))
  {
    tryCatch(search_info <- entrez_search(db="pubmed", term=search_term), error=function(e) NULL)
    
    if (is.null(search_info))
    {
      Sys.sleep(sleep)
      sleep <- sleep * attempt
      attempt <- attempt + 1
    }
  }
  
  search_summary <- NULL
  
  search_summary <- as.data.frame(paste0(as.character(Sys.time()), ", ",
                                         search_term, ", ",
                                         search_info$count))
  return(search_summary)
}

#' Estimate linear-model stats for a matrix of data using limma with empirical Bayes moderated t-stats for p-values
#'
#' @param mat: Nxp data matrix with N cell lines and p genes
#' @param vec: N vector of independent variables. Can be two-group labels as factors, bools, or can be numeric
#' @param covars: Optional Nxk matrix of covariates
#' @param weights: Optional N vector of precision weights for each data point
#' @param target_type: Name of the column variable in the data (default 'Gene')
#' @param limma_trend: Whether to fit an intensity trend with the empirical Bayes variance model
#'
#' @return: data frame of stats
#' @export
#'
#' @examples
#' CRISPR = load.from.taiga(data.name='avana-2-0-1-d98f',
#' data.version=1,
#' data.file='ceres_gene_effects',
#' transpose = T)
#' is_panc <- load.from.taiga(data.name = 'ccle-lines-lineages') %>% .[, 'pancreas']
#' ulines <- intersect(rownames(CRISPR), names(is_panc))
#' lim_res <- run_lm_stats_limma(CRISPR[ulines,], is_panc[ulines])
#' @export run_lm_stats_limma
run_lm_stats_limma <- function(mat, vec, covars = NULL, weights = NULL, target_type = 'Gene', limma_trend = FALSE) {
  require(limma)
  require(magrittr)
  require(tibble)
  require(plyr)
  require(dplyr)
  
  udata <- which(!is.na(vec))
  if (!is.numeric(vec)) {
    pred <- factor(vec[udata])
    stopifnot(length(levels(pred)) == 2) #only two group comparisons implemented so far
    n_out <- colSums(!is.na(mat[udata[pred == levels(pred)[1]],,drop=F]))
    n_in <- colSums(!is.na(mat[udata[pred == levels(pred)[2]],,drop=F]))
    min_samples <- pmin(n_out, n_in) %>% set_names(colnames(mat))
  } else {
    pred <- vec[udata]
    min_samples <- colSums(!is.na(mat[udata,]))
  }
  #there must be more than one unique value of the independent variable
  if (length(unique(pred)) <= 1) {
    return(NULL)
  }
  #if using covariates add them as additional predictors to the model
  if (!is.null(covars)) {
    if (!is.data.frame(covars)) {
      covars <- data.frame(covars)
    }
    combined <- covars[udata,, drop = FALSE]
    combined[['pred']] <- pred
    form <- as.formula(paste('~', paste0(colnames(combined), collapse = ' + ')))
    design <- model.matrix(form, combined)
    design <- design[, colSums(design) != 0, drop = FALSE]
  } else {
    design <- model.matrix(~pred)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata,])
    } else{
      weights <- weights[udata]
    }
  }
  fit <- limma::lmFit(t(mat[udata,]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- grep('pred', colnames(design), value = TRUE)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf)
  
  if (colnames(results)[1] == 'ID') {
    colnames(results)[1] <- target_type
  } else {
    results %<>% rownames_to_column(var = target_type)
  }
  results$min_samples <- min_samples[results[[target_type]]]
  
  two_to_one_sided <- function(two_sided_p, stat, test_dir) {
    #helper function for converting two-sided p-values to one-sided p-values
    one_sided_p <- two_sided_p / 2
    if (test_dir == 'right') {
      one_sided_p[stat < 0] <- 1 - one_sided_p[stat < 0]
    } else {
      one_sided_p[stat > 0] <- 1 - one_sided_p[stat > 0]
    }
    return(one_sided_p)
  }
  results %<>% set_colnames(revalue(colnames(.), c('logFC' = 'EffectSize', 'AveExpr' = 'Avg', 't' = 't_stat', 'B' = 'log_odds',
                                                   'P.Value' = 'p.value', 'adj.P.Val' = 'q.value', 'min_samples' = 'min_samples'))) %>% na.omit()
  results %<>% dplyr::mutate(p.left = two_to_one_sided(p.value, EffectSize, 'left'),
                             p.right = two_to_one_sided(p.value, EffectSize, 'right'),
                             q.left = p.adjust(p.left, method = 'BH'),
                             q.right = p.adjust(p.right, method = 'BH'))
  return(results)
}



#' Estimate linear-model stats for a matrix of data wrt a group of phenotype variables using limma with empirical Bayes moderated F-stats for p-values.
#' Note: Design matrix only models levels with at least 3 samples (the rest of samples are just modeled with the intercept term). Otherwise, limma hits errors
#'
#' @param mat: Nxp data matrix with N cell lines and p genes
#' @param phenos: Nxq dataframe/matrix of independent variables.
#' @param covars: Optional Nxk matrix of covariates
#' @param weights: Optional N vector of precision weights for each data point
#' @param target_type: Name of the column variable in the data (default 'Gene')
#' @param limma_trend: Whether to fit an intensity trend with the empirical Bayes variance model
#'
#' @return: data frame of stats
#' @export
#'
#' @examples
#' gene.effect <- load.from.taiga(data.name='avana-broad-18q1-ec8b', data.version=2, data.file='gene_effect')
#' sample.info <- load.from.taiga(data.name='avana-broad-18q1-ec8b', data.version=2, data.file='sample_info')
#'
#' covars <- sample.info %>%
#'   dplyr::select(CCLE_ID = cell_line, cas9_activity, cell_line_SSMD) %>%
#'   mutate(cas9_activity = as.numeric(cas9_activity)) %>%
#'   filter(!is.na(cas9_activity)) %>%
#'   column_to_rownames(var = 'CCLE_ID')
#' uvars <- sample.info %>%
#'   dplyr::select(CCLE_ID = cell_line, primary_tissue) %>%
#'   column_to_rownames(var = 'CCLE_ID')
#' res <- cdsr::run_lm_stats_limma_group(gene.effect, uvars, covars, weights = NULL, target_type = 'Gene')
#'
#' @export run_lm_stats_limma_group
run_lm_stats_limma_group <- function(mat, phenos, covars = NULL, weights = NULL, target_type = 'Gene', limma_trend = FALSE) {
  require(limma)
  require(magrittr)
  require(tibble)
  require(plyr)
  require(dplyr)
  
  udata <- rownames(mat) %>%
    intersect(rownames(phenos))
  if (!is.null(covars)) {
    udata %<>% intersect(rownames(covars))
  }
  
  form <- as.formula(paste('~', paste0(colnames(phenos), collapse = ' + ')))
  design <- model.matrix(form, data = phenos[udata,,drop=F])
  
  if (!is.null(covars)) {
    covars <- data.frame(covars)
    form <- as.formula(paste('~', paste0(colnames(covars), collapse = ' + ')))
    Cdesign <- model.matrix(form, data = covars[udata,,drop=F])
    Cdesign <- Cdesign[, setdiff(colnames(Cdesign), '(Intercept)'), drop = FALSE]
    stopifnot(length(intersect(colnames(Cdesign), colnames(design))) == 0)
    design %<>% cbind(Cdesign)
  }
  
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata,])
    } else{
      weights <- weights[udata]
    }
  }
  
  design <- design[, colSums(design) > 2, drop = FALSE] #require more than two samples from each group (otherwise returns errors)
  targ_coefs <- setdiff(colnames(design), '(Intercept)')
  fit <- limma::lmFit(t(mat[udata,]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- which(colnames(design) %in% targ_coefs)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf, sort.by = 'F', genelist = colnames(mat))
  results %<>% rownames_to_column(var = target_type)
  results %<>% set_colnames(revalue(colnames(.), c('AveExpr' = 'Avg',
                                                   'F' = 'F_stat',
                                                   'P.Value' = 'p.value',
                                                   'adj.P.Val' = 'q.value'))) %>%
    na.omit() %>%
    dplyr::select(-ProbeID)
  
  return(results)
}
