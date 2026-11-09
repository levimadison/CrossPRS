
#' Download Neale Lab UK Biobank GWAS Summary Statistics
#'
#' Downloads Neale Lab UK Biobank GWAS summary statistics for a specified
#' phenotype and one or more sex groups. The files are saved in the
#' destination directory in bgzip-compressed format, such as \code{.tsv.bgz}.
#'
#' Phenotype codes and corresponding summary-statistics download links are
#' listed in the Neale Lab UK Biobank GWAS Imputed v3 File Manifest Release
#' 20180731.
#'
#' @param pheno_code Character string. Phenotype code identifying the Neale
#'   Lab UK Biobank GWAS phenotype to download. For example,
#'   \code{"21001_irnt"} refers to inverse rank normalized body mass index
#'   (BMI), while \code{"21001_raw"} refers to untransformed BMI.
#'
#' @param sex_groups Character vector specifying the sex groups for which
#'   summary statistics should be downloaded. Supported values are:
#'   \describe{
#'     \item{\code{"both_sexes"}}{Results based on male and female
#'       participants.}
#'     \item{\code{"female"}}{Results based on female participants only.}
#'     \item{\code{"male"}}{Results based on male participants only.}
#'   }
#'   The default is \code{c("both_sexes", "female", "male")}. To download
#'   just one group, supply a single value, such as
#'   \code{sex_groups = "female"}.
#'
#' @param dest_dir Character string. Directory where the downloaded GWAS
#'   summary-statistics files will be saved. Defaults to the current working
#'   directory, \code{"./"}.
#'
#' @return Invisibly returns information about the downloaded files, if
#'   available. The downloaded GWAS summary-statistics files are saved in
#'   \code{dest_dir}.
#'
#' @details
#' ## Phenotype codes and download links
#'
#' The Neale Lab UK Biobank GWAS Imputed v3 File Manifest Release 20180731
#' contains a \code{README} sheet describing the available files and a
#' \code{Manifest 201807} sheet listing phenotype codes, phenotype
#' descriptions, and download links:
#'
#' \url{https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU}
#'
#' For example, the female-specific GWAS summary-statistics files for BMI are:
#' \itemize{
#'   \item Inverse rank normalized BMI, phenotype code
#'     \code{21001_irnt}:
#'     \url{https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_irnt.gwas.imputed_v3.female.tsv.bgz}
#'   \item Untransformed BMI, phenotype code \code{21001_raw}:
#'     \url{https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_raw.gwas.imputed_v3.female.tsv.bgz}
#' }
#'
#' These URLs illustrate the female-specific files. Use the manifest to find
#' the corresponding download link for other phenotypes and sex groups.
#'
#' The downloaded GWAS file can be supplied as \code{s_bgzfile} to
#' [build_prs_weights()] when constructing PRS SNP weights, or as
#' \code{t_bgzfile} to [test_prs_association()] when testing a PRS against a
#' target GWAS. The variant annotation file used by
#' [build_prs_weights()] is a separate file and is not described as a
#' download performed by this function.
#'
#' @examples
#' \dontrun{
#' ## Download female-specific GWAS summary statistics for inverse rank
#' ## normalized BMI into the "gwas_data" directory.
#' download_neale_gwas(
#'   pheno_code = "21001_irnt",
#'   sex_groups = "female",
#'   dest_dir = "gwas_data"
#' )
#'
#' ## Download summary statistics for all supported sex groups using the
#' ## default sex_groups setting.
#' download_neale_gwas(
#'   pheno_code = "21001_raw",
#'   dest_dir = "gwas_data"
#' )
#' }
#'
#' @seealso
#' [build_prs_weights()],
#' [test_prs_association()],
#' \url{https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU}
#'
#' @import curl
#' @export
download_neale_gwas<-function(pheno_code, sex_groups = c("both_sexes", "female", "male"), dest_dir="./") {
	# Base URL for Neale Lab Round 2 additive GWAS summary stats
	base_url <- "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/"
	
	for (sex in sex_groups) {
		
		# Biomarkers starting with '30' require .varorder. in the URL/filename
		is_varorder <- grepl("^30", pheno_code)
		
		suffix <- if (is_varorder) {
					paste0(".gwas.imputed_v3.", sex, ".varorder.tsv.bgz")
				} else {
					paste0(".gwas.imputed_v3.", sex, ".tsv.bgz")
				}
		
		file_name <- paste0(pheno_code, suffix)
		download_url <- paste0(base_url, file_name)
		dest_path <- file.path(dest_dir, file_name)
		
		# Skip if file already exists locally to avoid re-downloading
		if (file.exists(dest_path)) {
			cat(sprintf("Skipping (already downloaded): %s\n", file_name))
			next
		}
		
		cat(sprintf("Downloading: %s ...\n", file_name))
		
		tryCatch({
#					download.file(
#							url = download_url, 
#							destfile = dest_path, 
#							mode = "wb",          # Binary mode required for .bgz files
#							quiet = FALSE
#					)
					curl::curl_download(url= download_url, destfile= dest_path, quiet = FALSE)
					cat(sprintf("Successfully downloaded: %s\n\n", file_name))
				}, error = function(e) {
					cat(sprintf("ERROR downloading %s: %s\n\n", file_name, e$message))
				})
	}
}




