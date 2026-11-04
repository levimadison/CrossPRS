
#' Core Summary-Statistics PRS Association Method
#'
#' Implements the core summary-statistics polygenic risk score (PRS)
#' association method developed in this package. Given SNP weights, SNP
#' variances, SNP effect estimates, and their variances, the function
#' estimates the association between a PRS and a phenotype using GWAS summary
#' statistics.
#'
#' This is a low-level computational function intended for users who already
#' have the required summary statistics vectors. For analyses based on Neale
#' Lab UK Biobank GWAS summary statistics files, users should typically call
#' [test_prs_association()], which performs data input, SNP matching, allele
#' alignment, and then calls this function internally.
#'
#' @param w_i Numeric vector of SNP weights. Each element represents the
#'   weight assigned to a candidate independent SNP in the PRS.
#'
#' @param var_x Numeric vector of SNP variances. Each element represents the
#'   variance of the corresponding SNP.
#'
#' @param b_hat Numeric vector of SNP effect estimates. Each element
#'   represents the estimated effect (beta coefficient) of the corresponding
#'   SNP on the target phenotype.
#'
#' @param var_b_hat Numeric vector of variances of the SNP effect estimates.
#'   Each element corresponds to the variance of the associated value in
#'   `b_hat`.
#'
#' @return
#' A named numeric vector containing:
#' \describe{
#'   \item{effect_est}{Estimated effect of the PRS on the phenotype.}
#'   \item{effect_se}{Standard error of the estimated PRS effect.}
#'   \item{Z_a}{Z statistic for testing the null hypothesis that the PRS
#'   effect is zero.}
#'   \item{p_value}{Two-sided p-value corresponding to the Z statistic.}
#' }
#'
#' @details
#' This function implements the summary-statistics PRS association method
#' developed in this package. It estimates the effect of a PRS directly from
#' GWAS summary statistics without requiring individual-level genotype or
#' phenotype data.
#'
#' The function assumes that the input vectors correspond to the same set of
#' candidate independent SNPs and are ordered consistently. Specifically,
#' `w_i`, `var_x`, `b_hat`, and `var_b_hat` must have the same length, with
#' each element referring to the same SNP across all vectors.
#'
#' @examples
#' \dontrun{
#' ## Example vectors
#' w_i <- c(0.15, -0.08, 0.23)
#' var_x <- c(0.42, 0.39, 0.47)
#' b_hat <- c(0.031, -0.014, 0.052)
#' var_b_hat <- c(0.00012, 0.00010, 0.00015)
#'
#' prs_assoc_core(
#'   w_i = w_i,
#'   var_x = var_x,
#'   b_hat = b_hat,
#'   var_b_hat = var_b_hat
#' )
#' }
#'
#' @seealso
#' [build_prs_weights()],
#' [test_prs_association()]
#'
#' @export
prs_assoc_core<-function(w_i,var_x,  b_hat, var_b_hat) {
	#PRS effect estimator: alpha_hat
	effect_numerator <- sum(w_i * var_x * b_hat, na.rm = TRUE)
	effect_denominator <-sum((w_i^2) * (var_x))
	effect_est<-effect_numerator/effect_denominator
	
#PRS effect estimator variance: var(alpha_hat)
	effectVar_numerator<-sum((w_i^2) * (var_x^2) * var_b_hat, na.rm = TRUE)
	effectVar_denominator<-effect_denominator^2
	effectVar<-effectVar_numerator/effectVar_denominator
	effect_se<-sqrt(effectVar)
	
# Method 1: Calculate PRS statistic Z_a
#Z_a <- effect_est / effect_se
	
# PRS statistic Z Numerator: Sum_i [ w_i * Var(X_i) * b_hat_i ]
	z_numerator <- sum(w_i * var_x * b_hat, na.rm = TRUE)
	
# PRS statistic Z Denominator: sqrt( Sum_i [ w_i^2 * Var(X_i)^2 * Var(b_hat_i) ] )
	z_denominator <- sqrt(sum((w_i^2) * (var_x^2) * var_b_hat, na.rm = TRUE))
	
# Calculate PRS statistic Z_a
	Z_a <- z_numerator / z_denominator
	
# Compute two-tailed p-value from standard normal distribution
	p_value <- 2 * pnorm(-abs(Z_a))
	
	# Return named numeric vector of length 5
	res <- c(
			effect_est = as.numeric(effect_est),
			effect_se  = as.numeric(effect_se),
			Z_a        = as.numeric(Z_a),
			p_value    = as.numeric(p_value)
	)
	
	return(res)
}
