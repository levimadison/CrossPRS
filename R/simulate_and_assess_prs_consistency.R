#' Assess Consistency Between Two PRS Association Methods
#'
#' Simulate multiple independent datasets and assess the consistency
#' between an individual-level PRS association analysis and the
#' summary-statistic PRS association analysis implemented by
#' [prs_assoc_core()].
#'
#' For each simulated dataset, this function:
#' \enumerate{
#'   \item Simulates SNP genotypes and two outcomes using
#'         [simulate_snps_two_outcomes()].
#'   \item Performs SNP-level and PRS-level analyses using
#'         [analyze_snps_two_outcomes()].
#'   \item Extracts the results from analysis 3, which regresses `y2`
#'         on an individual-level PRS constructed using the SNP-level
#'         `z1` statistics.
#'   \item Extracts the results from analysis 4, which calls
#'         [prs_assoc_core()] using SNP-level summary statistics.
#'   \item Compares the two methods using correlation, bias, RMSE,
#'         maximum absolute difference, and agreement in statistical
#'         significance.
#' }
#'
#' The function also optionally generates scatter plots with identity
#' lines and Bland--Altman plots for the estimated effect and z
#' statistic.
#'
#' @param n_sim Integer. Number of independent datasets to simulate.
#'   Must be a positive integer. Default is `1000`.
#'
#' @param N Integer. Number of independent SNPs to simulate in each
#'   dataset. Passed to [simulate_snps_two_outcomes()].
#'
#' @param M Integer. Number of SNPs, among the `N` SNPs, that are
#'   associated with `y2`. Must satisfy `0 <= M <= N`.
#'
#' @param same_direction Logical. Passed to
#'   [simulate_snps_two_outcomes()]. If `TRUE`, the effects of the
#'   `M` SNPs on `y2` have the same direction as their effects on
#'   `y1`.
#'
#' @param n_ind Integer. Number of individuals in each simulated
#'   dataset. Passed to [simulate_snps_two_outcomes()]. Default is
#'   `1000`.
#'
#' @param effect_mean Numeric. Mean of the absolute SNP effect-size
#'   distribution. Passed to [simulate_snps_two_outcomes()]. Default
#'   is `0.10`.
#'
#' @param effect_sd Numeric. Standard deviation of the absolute SNP
#'   effect-size distribution. Passed to
#'   [simulate_snps_two_outcomes()]. Default is `0.02`.
#'
#' @param maf_min Numeric. Minimum minor allele frequency. Passed to
#'   [simulate_snps_two_outcomes()]. Must be greater than 0 and less
#'   than 0.50. Default is `0.01`.
#'
#' @param seed Optional integer used to initialize the random-number
#'   generator. If `NULL`, the current random-number state is used.
#'   Default is `NULL`.
#'
#' @param alpha Numeric significance level used to determine whether
#'   each method produces a statistically significant result.
#'   Default is `0.05`.
#'
#' @param make_plots Logical. If `TRUE`, create scatter/identity-line
#'   and Bland--Altman plots. Requires the `ggplot2` package.
#'   Default is `TRUE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{simulation_results}{
#'     A data frame containing the results from analysis 3 and
#'     analysis 4 for every simulated dataset. It contains the
#'     estimated effect, standard error, z statistic, p-value, and
#'     their pairwise differences.
#'   }
#'
#'   \item{agreement_summary}{
#'     A data frame summarizing agreement for the effect estimate,
#'     standard error, and z statistic. The summary includes the
#'     number of valid simulations, Pearson correlation, mean
#'     difference (bias), standard deviation of the difference, RMSE,
#'     mean absolute difference, and maximum absolute difference.
#'   }
#'
#'   \item{significance_agreement}{
#'     A contingency table comparing whether analysis 3 and analysis 4
#'     are statistically significant at the specified `alpha` level.
#'   }
#'
#'   \item{significance_concordance}{
#'     The proportion of simulations for which analysis 3 and analysis
#'     4 give the same significance conclusion.
#'   }
#'
#'   \item{plots}{
#'     A named list containing the generated `ggplot2` plots when
#'     `make_plots = TRUE`. The list contains identity-line plots for
#'     beta, SE, z, and p-value, and Bland--Altman plots for beta and z.
#'     If plots are not requested or `ggplot2` is unavailable, this
#'     component is `NULL`.
#' }
#' }
#'
#' @details
#' Analysis 3 constructs an individual-level PRS as
#'
#' \deqn{
#' PRS_i = \sum_{j=1}^{N} z_{1j}G_{ij},
#' }
#'
#' where `G` is the genotype dosage matrix and `z1` is the SNP-level
#' association statistic obtained from the analysis of `y1`.
#'
#' The third analysis then fits
#'
#' \deqn{
#' y_{2i} = \alpha + \beta_{PRS}PRS_i + \epsilon_i.
#' }
#'
#' The resulting effect estimate, standard error, z statistic, and
#' p-value are compared with the output of `prs_assoc_core()`.
#'
#' For analysis 4, `analyze_snps_two_outcomes()` calls
#' `prs_assoc_core()` using
#'
#' \itemize{
#'   \item `w_i = z1`,
#'   \item `var_x = 2 * maf * (1 - maf)`,
#'   \item `b_hat = beta2`, and
#'   \item `var_b_hat = se2^2`.
#' }
#'
#' The two methods are expected to give similar results when the
#' assumptions underlying the summary-statistic method hold.
#'
#' This function evaluates numerical agreement rather than relying
#' primarily on hypothesis tests. This is important because a paired
#' hypothesis test can detect very small differences when many
#' simulation replicates are used. Correlation, bias, RMSE, maximum
#' absolute difference, and identity-line plots provide more direct
#' measures of agreement.
#'
#' Small differences between the methods can occur because the
#' individual-level PRS analysis uses the observed genotype matrix,
#' whereas `prs_assoc_core()` uses the theoretical genotype variance
#' `2 * maf * (1 - maf)`.
#'
#' @examples
#' \dontrun{
#' # Run 1,000 simulation replicates.
#' consistency <- simulate_and_assess_prs_consistency(
#'   n_sim = 1000,
#'   N = 100,
#'   M = 30,
#'   same_direction = TRUE,
#'   n_ind = 5000,
#'   effect_mean = 0.10,
#'   effect_sd = 0.02,
#'   maf_min = 0.01,
#'   seed = 123
#' )
#'
#' # Quantitative agreement.
#' consistency$agreement_summary
#'
#' # Results from every simulation.
#' head(consistency$simulation_results)
#'
#' # Agreement of statistical significance at alpha = 0.05.
#' consistency$significance_agreement
#'
#' # Proportion of simulations with the same significance conclusion.
#' consistency$significance_concordance
#'
#' # Display plots.
#' consistency$plots$beta
#' consistency$plots$se
#' consistency$plots$z
#' consistency$plots$p_value
#'
#' # Bland--Altman plots.
#' consistency$plots$beta_bland_altman
#' consistency$plots$z_bland_altman
#' }
#'
#' @export
simulate_and_assess_prs_consistency <- function(
		n_sim = 1000,
		N,
		M,
		same_direction,
		n_ind = 1000,
		effect_mean = 0.10,
		effect_sd = 0.02,
		maf_min = 0.01,
		seed = NULL,
		alpha = 0.05,
		make_plots = TRUE) {
	
	# ============================================================
	# Input validation
	# ============================================================
	
	if (!is.numeric(n_sim) ||
			length(n_sim) != 1L ||
			!is.finite(n_sim) ||
			n_sim <= 0 ||
			n_sim != as.integer(n_sim)) {
		stop("n_sim must be a positive integer.")
	}
	
	if (!is.numeric(N) ||
			length(N) != 1L ||
			!is.finite(N) ||
			N <= 0 ||
			N != as.integer(N)) {
		stop("N must be a positive integer.")
	}
	
	if (!is.numeric(M) ||
			length(M) != 1L ||
			!is.finite(M) ||
			M < 0 ||
			M > N ||
			M != as.integer(M)) {
		stop("M must be an integer between 0 and N.")
	}
	
	if (!is.logical(same_direction) ||
			length(same_direction) != 1L ||
			is.na(same_direction)) {
		stop("same_direction must be TRUE or FALSE.")
	}
	
	if (!is.numeric(n_ind) ||
			length(n_ind) != 1L ||
			!is.finite(n_ind) ||
			n_ind <= 0 ||
			n_ind != as.integer(n_ind)) {
		stop("n_ind must be a positive integer.")
	}
	
	if (!is.numeric(effect_mean) ||
			length(effect_mean) != 1L ||
			!is.finite(effect_mean) ||
			effect_mean <= 0) {
		stop("effect_mean must be a positive number.")
	}
	
	if (!is.numeric(effect_sd) ||
			length(effect_sd) != 1L ||
			!is.finite(effect_sd) ||
			effect_sd < 0) {
		stop("effect_sd must be non-negative.")
	}
	
	if (!is.numeric(maf_min) ||
			length(maf_min) != 1L ||
			!is.finite(maf_min) ||
			maf_min <= 0 ||
			maf_min >= 0.50) {
		stop("maf_min must be greater than 0 and less than 0.50.")
	}
	
	if (!is.null(seed)) {
		if (!is.numeric(seed) ||
				length(seed) != 1L ||
				!is.finite(seed)) {
			stop("seed must be NULL or a single finite number.")
		}
	}
	
	if (!is.numeric(alpha) ||
			length(alpha) != 1L ||
			!is.finite(alpha) ||
			alpha <= 0 ||
			alpha >= 1) {
		stop("alpha must be between 0 and 1.")
	}
	
	if (!is.logical(make_plots) ||
			length(make_plots) != 1L ||
			is.na(make_plots)) {
		stop("make_plots must be TRUE or FALSE.")
	}
	
	# Convert integer-like parameters to integer.
	n_sim <- as.integer(n_sim)
	N <- as.integer(N)
	M <- as.integer(M)
	n_ind <- as.integer(n_ind)
	
	# Set random-number seed.
	if (!is.null(seed)) {
		set.seed(seed)
	}
	
	# ============================================================
	# Storage for simulation results
	# ============================================================
	
	simulation_results <- data.frame(
			simulation = seq_len(n_sim),
			
			# Analysis 3: individual-level PRS
			beta_prs = rep(NA_real_, n_sim),
			se_prs = rep(NA_real_, n_sim),
			z_prs = rep(NA_real_, n_sim),
			p_prs = rep(NA_real_, n_sim),
			
			# Analysis 4: prs_assoc_core()
			effect_est = rep(NA_real_, n_sim),
			effect_se = rep(NA_real_, n_sim),
			Z_a = rep(NA_real_, n_sim),
			p_core = rep(NA_real_, n_sim),
			
			# Differences: analysis 3 - analysis 4
			diff_beta = rep(NA_real_, n_sim),
			diff_se = rep(NA_real_, n_sim),
			diff_z = rep(NA_real_, n_sim)
	)
	
	# ============================================================
	# Run simulations
	# ============================================================
	
	for (i in seq_len(n_sim)) {
		
		# ----------------------------------------------------------
		# Simulate one dataset
		# ----------------------------------------------------------
		
		sim_data <- simulate_snps_two_outcomes(
				N = N,
				M = M,
				same_direction = same_direction,
				n_ind = n_ind,
				effect_mean = effect_mean,
				effect_sd = effect_sd,
				maf_min = maf_min
		)
		
		# ----------------------------------------------------------
		# Analyze simulated dataset
		# ----------------------------------------------------------
		
		analysis <- analyze_snps_two_outcomes(
				data = sim_data
		)
		
		# ----------------------------------------------------------
		# Extract analysis 3:
		# Individual-level PRS analysis
		# ----------------------------------------------------------
		
		simulation_results$beta_prs[i] <-
				unname(
						analysis$prs$association["beta"]
				)
		
		simulation_results$se_prs[i] <-
				unname(
						analysis$prs$association["se"]
				)
		
		simulation_results$z_prs[i] <-
				unname(
						analysis$prs$association["z"]
				)
		
		simulation_results$p_prs[i] <-
				unname(
						analysis$prs$association["p_value"]
				)
		
		# ----------------------------------------------------------
		# Extract analysis 4:
		# prs_assoc_core()
		# ----------------------------------------------------------
		
		simulation_results$effect_est[i] <-
				unname(
						analysis$prs_assoc_core["effect_est"]
				)
		
		simulation_results$effect_se[i] <-
				unname(
						analysis$prs_assoc_core["effect_se"]
				)
		
		simulation_results$Z_a[i] <-
				unname(
						analysis$prs_assoc_core["Z_a"]
				)
		
		simulation_results$p_core[i] <-
				unname(
						analysis$prs_assoc_core["p_value"]
				)
		
		# ----------------------------------------------------------
		# Calculate differences
		# ----------------------------------------------------------
		
		simulation_results$diff_beta[i] <-
				simulation_results$beta_prs[i] -
				simulation_results$effect_est[i]
		
		simulation_results$diff_se[i] <-
				simulation_results$se_prs[i] -
				simulation_results$effect_se[i]
		
		simulation_results$diff_z[i] <-
				simulation_results$z_prs[i] -
				simulation_results$Z_a[i]
	}
	
	# ============================================================
	# Helper function for agreement statistics
	# ============================================================
	
	calculate_agreement <- function(x, y) {
		
		keep <- is.finite(x) & is.finite(y)
		
		x <- x[keep]
		y <- y[keep]
		
		if (length(x) == 0L) {
			return(
					c(
							n = 0,
							correlation = NA_real_,
							mean_difference = NA_real_,
							sd_difference = NA_real_,
							RMSE = NA_real_,
							mean_absolute_difference = NA_real_,
							max_absolute_difference = NA_real_
					)
			)
		}
		
		difference <- x - y
		
		# Correlation is undefined when one or both vectors
		# have zero variance.
		correlation <- if (
						length(x) > 1L &&
						sd(x) > 0 &&
						sd(y) > 0
						) {
					cor(x, y)
				} else {
					NA_real_
				}
		
		c(
				n = length(x),
				
				correlation = correlation,
				
				mean_difference = mean(
						difference
				),
				
				sd_difference = if (length(difference) > 1L) {
							sd(difference)
						} else {
							NA_real_
						},
				
				RMSE = sqrt(
						mean(difference^2)
				),
				
				mean_absolute_difference = mean(
						abs(difference)
				),
				
				max_absolute_difference = max(
						abs(difference)
				)
		)
	}
	
	# ============================================================
	# Calculate agreement for beta, SE, and Z
	# ============================================================
	
	beta_agreement <- calculate_agreement(
			x = simulation_results$beta_prs,
			y = simulation_results$effect_est
	)
	
	se_agreement <- calculate_agreement(
			x = simulation_results$se_prs,
			y = simulation_results$effect_se
	)
	
	z_agreement <- calculate_agreement(
			x = simulation_results$z_prs,
			y = simulation_results$Z_a
	)
	
	agreement_summary <- rbind(
			beta = beta_agreement,
			se = se_agreement,
			z = z_agreement
	)
	
	agreement_summary <- as.data.frame(
			agreement_summary,
			stringsAsFactors = FALSE
	)
	
	agreement_summary$measure <-
			rownames(agreement_summary)
	
	rownames(agreement_summary) <- NULL
	
	agreement_summary <- agreement_summary[
			,
			c(
					"measure",
					"n",
					"correlation",
					"mean_difference",
					"sd_difference",
					"RMSE",
					"mean_absolute_difference",
					"max_absolute_difference"
			)
	]
	
	# ============================================================
	# Compare statistical significance
	# ============================================================
	
	valid_p <- is.finite(
					simulation_results$p_prs
			) &
			is.finite(
					simulation_results$p_core
			)
	
	sig_prs <- simulation_results$p_prs[
			valid_p
	] < alpha
	
	sig_core <- simulation_results$p_core[
			valid_p
	] < alpha
	
	significance_agreement <- table(
			Individual_PRS = sig_prs,
			PRS_Core = sig_core
	)
	
	significance_concordance <- if (
					length(sig_prs) > 0
					) {
				mean(
						sig_prs == sig_core
				)
			} else {
				NA_real_
			}
	
	# ============================================================
	# Create plots
	# ============================================================
	
	plots <- NULL
	
	if (make_plots) {
		
		if (!requireNamespace(
				"ggplot2",
				quietly = TRUE
		)) {
			
			warning(
					"ggplot2 is not installed. ",
					"No plots were created. ",
					"Install ggplot2 or use make_plots = FALSE."
			)
			
		} else {
			
			# --------------------------------------------------------
			# 1. Beta identity-line plot
			# --------------------------------------------------------
			
			plots$beta <- ggplot2::ggplot(
							simulation_results,
							ggplot2::aes(
									x = beta_prs,
									y = effect_est
							)
					) +
					ggplot2::geom_point(
							alpha = 0.6
					) +
					ggplot2::geom_abline(
							intercept = 0,
							slope = 1,
							color = "red",
							linewidth = 1
					) +
					ggplot2::labs(
							x = "Individual-level PRS beta",
							y = "prs_assoc_core effect estimate",
							title = "PRS Effect Estimate Agreement"
					) +
					ggplot2::theme_bw()
			
			# --------------------------------------------------------
			# 2. SE identity-line plot
			# --------------------------------------------------------
			
			plots$se <- ggplot2::ggplot(
							simulation_results,
							ggplot2::aes(
									x = se_prs,
									y = effect_se
							)
					) +
					ggplot2::geom_point(
							alpha = 0.6
					) +
					ggplot2::geom_abline(
							intercept = 0,
							slope = 1,
							color = "red",
							linewidth = 1
					) +
					ggplot2::labs(
							x = "Individual-level PRS SE",
							y = "prs_assoc_core SE",
							title = "PRS Standard Error Agreement"
					) +
					ggplot2::theme_bw()
			
			# --------------------------------------------------------
			# 3. Z identity-line plot
			# --------------------------------------------------------
			
			plots$z <- ggplot2::ggplot(
							simulation_results,
							ggplot2::aes(
									x = z_prs,
									y = Z_a
							)
					) +
					ggplot2::geom_point(
							alpha = 0.6
					) +
					ggplot2::geom_abline(
							intercept = 0,
							slope = 1,
							color = "red",
							linewidth = 1
					) +
					ggplot2::labs(
							x = "Individual-level PRS Z statistic",
							y = "prs_assoc_core Z statistic",
							title = "PRS Z Statistic Agreement"
					) +
					ggplot2::theme_bw()
			
			# --------------------------------------------------------
			# 4. P-value identity-line plot
			#
			# Use -log10(P) because raw p-values can span many orders
			# of magnitude.
			# --------------------------------------------------------
			
			p_plot_data <- simulation_results
			
			# Avoid Inf when p = 0 due to numerical underflow.
			p_plot_data$p_prs_plot <- pmax(
					p_plot_data$p_prs,
					.Machine$double.xmin
			)
			
			p_plot_data$p_core_plot <- pmax(
					p_plot_data$p_core,
					.Machine$double.xmin
			)
			
			plots$p_value <- ggplot2::ggplot(
							p_plot_data,
							ggplot2::aes(
									x = -log10(p_prs_plot),
									y = -log10(p_core_plot)
							)
					) +
					ggplot2::geom_point(
							alpha = 0.6
					) +
					ggplot2::geom_abline(
							intercept = 0,
							slope = 1,
							color = "red",
							linewidth = 1
					) +
					ggplot2::labs(
							x = "-log10(P), individual-level PRS",
							y = "-log10(P), prs_assoc_core",
							title = "PRS P-value Agreement"
					) +
					ggplot2::theme_bw()
			
			# --------------------------------------------------------
			# 5. Bland--Altman plot for beta
			# --------------------------------------------------------
			
			beta_mean <- (
						simulation_results$beta_prs +
						simulation_results$effect_est
						) / 2
			
			beta_diff <- (
						simulation_results$beta_prs -
						simulation_results$effect_est
						)
			
			beta_bias <- mean(
					beta_diff,
					na.rm = TRUE
			)
			
			beta_sd <- sd(
					beta_diff,
					na.rm = TRUE
			)
			
			beta_upper <- beta_bias +
					1.96 * beta_sd
			
			beta_lower <- beta_bias -
					1.96 * beta_sd
			
			beta_ba_data <- data.frame(
					mean = beta_mean,
					difference = beta_diff
			)
			
			plots$beta_bland_altman <-
					ggplot2::ggplot(
							beta_ba_data,
							ggplot2::aes(
									x = mean,
									y = difference
							)
					) +
					ggplot2::geom_point(
							alpha = 0.6
					) +
					ggplot2::geom_hline(
							yintercept = beta_bias,
							color = "red",
							linewidth = 1
					) +
					ggplot2::geom_hline(
							yintercept = c(
									beta_lower,
									beta_upper
							),
							linetype = "dashed"
					) +
					ggplot2::labs(
							x = "Mean of two beta estimates",
							y = "Difference: Analysis 3 - Analysis 4",
							title = "Bland--Altman Plot for PRS Effect"
					) +
					ggplot2::theme_bw()
			
			# --------------------------------------------------------
			# 6. Bland--Altman plot for SE
			# --------------------------------------------------------
			
			se_mean <- (
						simulation_results$se_prs +
						simulation_results$effect_se
						) / 2
			
			se_diff <- (
						simulation_results$se_prs -
						simulation_results$effect_se
						)
			
			se_bias <- mean(
					se_diff,
					na.rm = TRUE
			)
			
			se_sd <- sd(
					se_diff,
					na.rm = TRUE
			)
			
			se_upper <- se_bias +
					1.96 * se_sd
			
			se_lower <- se_bias -
					1.96 * se_sd
			
			se_ba_data <- data.frame(
					mean = se_mean,
					difference = se_diff
			)
			
			plots$se_bland_altman <-
					ggplot2::ggplot(
							se_ba_data,
							ggplot2::aes(
									x = mean,
									y = difference
							)
					) +
					ggplot2::geom_point(
							alpha = 0.6
					) +
					ggplot2::geom_hline(
							yintercept = se_bias,
							color = "red",
							linewidth = 1
					) +
					ggplot2::geom_hline(
							yintercept = c(
									se_lower,
									se_upper
							),
							linetype = "dashed"
					) +
					ggplot2::labs(
							x = "Mean of two SE estimates",
							y = "Difference: Analysis 3 - Analysis 4",
							title = "Bland--Altman Plot for PRS SE"
					) +
					ggplot2::theme_bw()
			
			# --------------------------------------------------------
			# 7. Bland--Altman plot for Z
			# --------------------------------------------------------
			
			z_mean <- (
						simulation_results$z_prs +
						simulation_results$Z_a
						) / 2
			
			z_diff <- (
						simulation_results$z_prs -
						simulation_results$Z_a
						)
			
			z_bias <- mean(
					z_diff,
					na.rm = TRUE
			)
			
			z_sd <- sd(
					z_diff,
					na.rm = TRUE
			)
			
			z_upper <- z_bias +
					1.96 * z_sd
			
			z_lower <- z_bias -
					1.96 * z_sd
			
			z_ba_data <- data.frame(
					mean = z_mean,
					difference = z_diff
			)
			
			plots$z_bland_altman <-
					ggplot2::ggplot(
							z_ba_data,
							ggplot2::aes(
									x = mean,
									y = difference
							)
					) +
					ggplot2::geom_point(
							alpha = 0.6
					) +
					ggplot2::geom_hline(
							yintercept = z_bias,
							color = "red",
							linewidth = 1
					) +
					ggplot2::geom_hline(
							yintercept = c(
									z_lower,
									z_upper
							),
							linetype = "dashed"
					) +
					ggplot2::labs(
							x = "Mean of two Z statistics",
							y = "Difference: Analysis 3 - Analysis 4",
							title = "Bland--Altman Plot for PRS Z Statistic"
					) +
					ggplot2::theme_bw()
		}
	}
	
	# ============================================================
	# Return results
	# ============================================================
	
	return(
			list(
					simulation_results = simulation_results,
					
					agreement_summary = agreement_summary,
					
					significance_agreement = significance_agreement,
					
					significance_concordance = significance_concordance,
					
					plots = plots
			)
	)
}
