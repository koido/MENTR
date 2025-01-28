#! /opt/homebrew/bin/Rscript

#--------------
#
# Calculate predictive accuracy
#
#--------------

#require(methods)
require(argparse)
parser <- ArgumentParser()
parser$add_argument("--input_file", required = T,
                    help = "Path for input file. Header is required. [*_type.test.expr_pred.txt.gz]")
parser$add_argument("--col_types", required = T,
                    help = "Input file's col_types for readr:read_delim function")
parser$add_argument("--out_dir_cmn", required = T,
                    help = "Path for output directory.")
parser$add_argument("--out_name", required = T,
                    help = "Output file name.")
parser$add_argument("--header_prediction", default = "pred",
                    help = "Column name for predicted values.")
parser$add_argument("--header_expression", default = "auto",
                    help = "Column name for expressino levels. Default: expr for AUROC analysis and log_expr for Spearman's rho analysis")
parser$add_argument("--threshold_binary", default = 0,
                    help = "Threshold for binalizing. Default: >0 is expressed.")
parser$add_argument("--auroc", action = 'store_true',
                    help = "If true, AUROC is evaluated. Otherwise, Spearman's rho is evaluated.")

args <- parser$parse_args()

# library for analysis
require(tidyverse)

# summarize function
## AUROC
s_fn_auroc <- function(.df, stats){
  
  .n <- .df %>% nrow
  .n_case <- .df %>% dplyr::filter(binary == 1) %>% nrow
  
  # exclude metric with #N_case < 10
  if(.n_case < 10){
    .stats <- data.frame(roc = NA, roc_95CI_L = NA, roc_95CI_H = NA,
                         n = .n, n_case = .n_case,
                         remarks = 'NA due to #case<10')
    stats <- rbind(stats, .stats)
    return(stats)   
  }
  
  # exception
  if(.n == .n_case){
    .stats <- data.frame(roc = NA, roc_95CI_L = NA, roc_95CI_H = NA,
                         n = .n, n_case = .n_case,
                         remarks = 'All records are the case')
    stats <- rbind(stats, .stats)
    return(stats)  
  }
  
  # ROC evaluation
  ROC <- roc(response = .df %>% dplyr::select(binary) %>% unlist %>% as.vector,
             predictor = .df %>% dplyr::select(pred) %>% unlist %>% as.vector,
             direction = "<",
             algorithm = 3)
  ROC_CI <- ci(ROC, method = 'delong')
  
  .roc <- ROC$auc %>% as.vector
  .roc_95CI_L <- ROC_CI[1]
  .roc_95CI_H <- ROC_CI[3]
  .stats <- data.frame(roc = .roc, roc_95CI_L = .roc_95CI_L, roc_95CI_H = .roc_95CI_H,
                       n = .n, n_case = .n_case, remarks = NA)
  stats <- rbind(stats, .stats)
  return(stats)    
}
## Spearman's rho
s_fn_spearman <- function(.df, stats){
  
  .n <- .df %>% nrow
  
  # exclude metric with #records < 10
  if(nrow(.df) < 10){
    .stats <- data.frame(r = NA, r_95CI_L = NA, r_95CI_H = NA,
                         p = NA, t = NA, df = NA, n = .n,
                         remarks = 'NA due to #case<10')
    stats <- rbind(stats, .stats)
    return(stats)   
  }
  
  # cor.test(rank(a), rank(b), method = "pearson") is equivalent with spearman, which can calculate 95% CI.
  # tie.method = average (default)
  .cor_res <- cor.test(x = rank(.df$pred), y = rank(.df$log_expr), method = "pearson")
  .r <- .cor_res$estimate
  .r_95CI_L <- .cor_res$conf.int[1]
  .r_95CI_H <- .cor_res$conf.int[2]
  .p = .cor_res$p.value
  .t <- .cor_res$statistic
  .free <- .cor_res$parameter
  .n <- .df %>% nrow
  .stats <- data.frame(r = .r, r_95CI_L = .r_95CI_L, r_95CI_H = .r_95CI_H,
                       p = .p, t = .t, df = .free, n = .n, remarks = NA)
  stats <- rbind(stats, .stats)
  return(stats)    
}

# read prediction file
pred <- read_delim(file = args$input_file, delim = "\t", col_types = args$col_types)
colnames(pred)[colnames(pred) == args$header_prediction] <- "pred"

# calculate
stats <- data.frame()
if(args$auroc){
  require(pROC)
  if(args$header_expression == "auto"){
    # remove missing
    pred <- pred %>% dplyr::filter(!is.na(expr))
    pred <- dplyr::mutate(pred, binary = ifelse(expr > args$threshold_binary, 1, 0))
  }else{
    # remove missing
    eval(parse(text = sprintf("pred <- dplyr::filter(pred, !is.na(%s))", args$header_expression)))
    eval(parse(text = sprintf("pred <- dplyr::mutate(pred, binary = ifelse(%s > %s, 1, 0))", args$header_expression, args$threshold_binary)))
  }
  stats <- s_fn_auroc(pred, stats)
}else{
  if(args$header_expression != "auto"){
    colnames(pred)[colnames(pred) == args$header_expression] <- "log_expr"
  }
  stats <- s_fn_spearman(pred, stats)
}

# save
out_file = file.path(args$out_dir_cmn, "each", args$out_name)
dir.create(file.path(args$out_dir_cmn, "each"), showWarnings = F, recursive = T)
write.table(stats, file = out_file, sep = "\t", quote = FALSE, row.names = F)
