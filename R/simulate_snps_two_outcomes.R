#' Simulate SNPs with Two Outcomes
#'
#' Simulate independent SNP genotypes and two quantitative outcomes (`y1`
#' and `y2`) under a specified pattern of SNP associations.
#'
#' All `N` simulated SNPs have nonzero effects on `y1`. Exactly `M` of
#' these SNPs are randomly selected to have nonzero effects on `y2`.
#' When `same_direction = TRUE`, the effects of the `M` SNPs on `y2`
#' have the same direction (sign) as their effects on `y1`. When
#' `same_direction = FALSE`, the direction of each effect on `y2` is
#' independently sampled.
#'
#' Minor allele frequencies (MAFs) are independently sampled from a
#' uniform distribution between `maf_min` and 0.50. Genotypes are
#' generated under Hardy-Weinberg equilibrium using a binomial model,
#' with dosage equal to the number of minor/alternate alleles.
#'
#' The outcomes are generated according to
#' \deqn{
#' y_1 = G\beta_1 + \epsilon_1
#' }
#' and
#' \deqn{
#' y_2 = G\beta_2 + \epsilon_2,
#' }
#' where `G` is the genotype matrix and the residuals are normally
#' distributed. The resulting outcomes are standardized to mean zero
#' and standard deviation one.
#'
#' @param N Integer. Number of independent SNPs to simulate. Must be
#'   greater than zero.
#' @param M Integer. Number of SNPs, among the `N` SNPs, that are
#'   associated with `y2`. Must satisfy `0 <= M <= N`.
#' @param same_direction Logical. If `TRUE`, the `M` SNP effects on
#'   `y2` have the same direction as their effects on `y1`. If `FALSE`,
#'   their directions are independently determined.
#' @param n_ind Integer. Number of individuals to simulate.
#'   Default is 1000.
#' @param effect_mean Numeric. Mean of the absolute SNP effect-size
#'   distribution. Default is 0.10.
#' @param effect_sd Numeric. Standard deviation of the absolute
#'   SNP effect-size distribution. Default is 0.02.
#' @param maf_min Numeric. Minimum minor allele frequency. Must be
#'   greater than 0 and less than 0.50. Default is 0.01.
#' @param seed Optional integer. Random-number seed for reproducibility.
#'   Default is `NULL`.
#'
#' @return A list containing:
#' \describe{
#'   \item{G}{An `n_ind` by `N` genotype matrix. Each column represents
#'   a SNP and contains dosages 0, 1, or 2.}
#'   \item{y1}{Standardized quantitative outcome associated with all
#'   `N` SNPs.}
#'   \item{y2}{Standardized quantitative outcome associated with exactly
#'   `M` SNPs.}
#'   \item{maf}{Vector of length `N` containing the minor allele
#'   frequency for each SNP.}
#'   \item{beta_y1}{Vector of length `N` containing SNP effects on `y1`.}
#'   \item{beta_y2}{Vector of length `N` containing SNP effects on
#'   `y2`. Zero indicates no association with `y2`.}
#'   \item{associated_y1}{Indices of SNPs associated with `y1`. This is
#'   always `1:N`.}
#'   \item{associated_y2}{Indices of the `M` SNPs associated with `y2`.}
#' }
#'
#' @details
#' The alternate allele is defined as the minor allele. For SNP `j`,
#' the genotype dosage is generated as
#' \deqn{G_j \sim Binomial(2, MAF_j).}
#'
#' The absolute effect sizes are sampled from a normal distribution
#' and their signs are independently assigned. For `y1`, all SNPs
#' receive nonzero effects. For `y2`, exactly `M` SNPs receive effects.
#'
#' If `same_direction = TRUE`, the sign of each selected `y2` effect
#' is copied from the corresponding `y1` effect. If `same_direction =
#' FALSE`, the sign is independently generated.
#'
#' Note that standardizing the simulated outcomes gives mean zero and
#' variance one, but because the outcomes contain genotype-dependent
#' genetic components, they are not mathematically guaranteed to have
#' an exactly standard normal marginal distribution. The residual
#' component is normally distributed.
#'
#' @examples
#' # Simulate 100 SNPs, with 30 associated with y2.
#' sim <- simulate_snps_two_outcomes(
#'   N = 100,
#'   M = 30,
#'   same_direction = TRUE,
#'   n_ind = 5000,
#'   seed = 123
#' )
#'
#' # Genotype matrix
#' dim(sim$G)
#'
#' # Number of SNPs associated with each outcome
#' length(sim$associated_y1)
#' length(sim$associated_y2)
#'
#' # Check that exactly 30 y2 effects are nonzero
#' sum(sim$beta_y2 != 0)
#'
#' # Check direction concordance for shared SNPs
#' all(
#'   sign(sim$beta_y1[sim$associated_y2]) ==
#'     sign(sim$beta_y2[sim$associated_y2])
#' )
#'
#' @export
simulate_snps_two_outcomes <- function(
		N,
		M,
		same_direction,
		n_ind = 1000,
		effect_mean = 0.10,
		effect_sd = 0.02,
		maf_min = 0.01,
		seed = NULL) {
	
	# ------------------------------------------------------------
	# Input validation
	# ------------------------------------------------------------
	
	if (!is.numeric(N) || length(N) != 1L ||
			!is.finite(N) || N <= 0 || N != as.integer(N)) {
		stop("N must be a positive integer.")
	}
	
	if (!is.numeric(M) || length(M) != 1L ||
			!is.finite(M) || M < 0 || M > N || M != as.integer(M)) {
		stop("M must be an integer between 0 and N.")
	}
	
	if (!is.logical(same_direction) ||
			length(same_direction) != 1L ||
			is.na(same_direction)) {
		stop("same_direction must be TRUE or FALSE.")
	}
	
	if (!is.numeric(n_ind) || length(n_ind) != 1L ||
			!is.finite(n_ind) || n_ind <= 0 ||
			n_ind != as.integer(n_ind)) {
		stop("n_ind must be a positive integer.")
	}
	
	if (!is.numeric(effect_mean) || length(effect_mean) != 1L ||
			!is.finite(effect_mean) || effect_mean <= 0) {
		stop("effect_mean must be a positive number.")
	}
	
	if (!is.numeric(effect_sd) || length(effect_sd) != 1L ||
			!is.finite(effect_sd) || effect_sd < 0) {
		stop("effect_sd must be non-negative.")
	}
	
	if (!is.numeric(maf_min) || length(maf_min) != 1L ||
			!is.finite(maf_min) || maf_min <= 0 || maf_min >= 0.50) {
		stop("maf_min must be between 0 and 0.50.")
	}
	
	if (!is.null(seed)) {
		if (!is.numeric(seed) || length(seed) != 1L ||
				!is.finite(seed)) {
			stop("seed must be NULL or a single finite number.")
		}
		set.seed(seed)
	}
	
	N <- as.integer(N)
	M <- as.integer(M)
	n_ind <- as.integer(n_ind)
	
	# ------------------------------------------------------------
	# 1. Simulate minor allele frequencies
	# ------------------------------------------------------------
	
	maf <- runif(
			N,
			min = maf_min,
			max = 0.50
	)
	
	# ------------------------------------------------------------
	# 2. Simulate genotype dosages
	#
	# G_ij = number of minor/alternate alleles
	# ------------------------------------------------------------
	
	G <- matrix(
			nrow = n_ind,
			ncol = N
	)
	
	for (j in seq_len(N)) {
		G[, j] <- rbinom(
				n = n_ind,
				size = 2,
				prob = maf[j]
		)
	}
	
	colnames(G) <- paste0("SNP", seq_len(N))
	
	# ------------------------------------------------------------
	# 3. Generate effects on y1
	#
	# All N SNPs are associated with y1.
	# ------------------------------------------------------------
	
	beta_y1_abs <- abs(
			rnorm(
					N,
					mean = effect_mean,
					sd = effect_sd
			)
	)
	
	# Make sure effects are nonzero.
	beta_y1_abs[beta_y1_abs == 0] <- effect_mean
	
	beta_y1 <- beta_y1_abs * sample(
			c(-1, 1),
			size = N,
			replace = TRUE
	)
	
	
	# ------------------------------------------------------------
	# 4. Select M SNPs associated with y2
	# ------------------------------------------------------------
	
	associated_y2 <- if (M > 0) {
				sample(
						seq_len(N),
						size = M,
						replace = FALSE
				)
			} else {
				integer(0)
			}
	
	beta_y2 <- numeric(N)
	
	if (M > 0) {
		
		# Generate absolute effect sizes for y2
		beta_y2_abs <- abs(
				rnorm(
						M,
						mean = effect_mean,
						sd = effect_sd
				)
		)
		
		beta_y2_abs[beta_y2_abs == 0] <- effect_mean
		
		if (same_direction) {
			
			# Same sign as the corresponding y1 effect
			beta_y2[associated_y2] <-
					beta_y2_abs *
					sign(beta_y1[associated_y2])
			
		} else {
			
			# Independent direction
			beta_y2[associated_y2] <-
					beta_y2_abs *
					sample(
							c(-1, 1),
							size = M,
							replace = TRUE
					)
		}
	}
	
	# ------------------------------------------------------------
	# 5. Generate genetic components
	# ------------------------------------------------------------
	
	genetic_y1 <- as.vector(G %*% beta_y1)
	genetic_y2 <- as.vector(G %*% beta_y2)
	
	# ------------------------------------------------------------
	# 6. Generate normally distributed residuals
	# ------------------------------------------------------------
	
	epsilon_y1 <- rnorm(n_ind)
	epsilon_y2 <- rnorm(n_ind)
	
	y1 <- genetic_y1 + epsilon_y1
	y2 <- genetic_y2 + epsilon_y2
	
	# ------------------------------------------------------------
	# 7. Standardize outcomes
	# ------------------------------------------------------------
	
	y1 <- as.numeric(scale(y1))
	y2 <- as.numeric(scale(y2))
	
	# ------------------------------------------------------------
	# 8. Return simulation results
	# ------------------------------------------------------------
	
	return(list(
					G = G,
					y1 = y1,
					y2 = y2,
					maf = maf,
					beta_y1 = beta_y1,
					beta_y2 = beta_y2,
					associated_y1 = seq_len(N),
					associated_y2 = associated_y2
			))
}

