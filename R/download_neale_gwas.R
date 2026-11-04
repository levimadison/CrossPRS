
#' Download Neale Lab UK Biobank GWAS Summary Statistics
#'
#' Downloads Neale Lab UK Biobank GWAS summary statistics for a specified
#' phenotype and one or more sex groups. The summary statistics files are
#' downloaded in bgzip-compressed format (e.g., `.bgz` or `.bgz.gz`) and
#' saved to the specified destination directory.
#'
#' The function supports downloading GWAS summary statistics generated from
#' the UK Biobank GWAS Imputed v3 resource. By default, summary statistics
#' are downloaded for all available sex groups: both sexes, females only,
#' and males only. Users can restrict downloads to specific sex groups using
#' the `sex_groups` argument.
#'
#' @param pheno_code Character string. Phenotype code identifying the Neale
#'   Lab UK Biobank GWAS phenotype to download. For example,
#'   `"100011_raw"` identifies a specific phenotype in the UK Biobank GWAS
#'   resource.
#'
#' @param sex_groups Character vector specifying the sex groups for which
#'   GWAS summary statistics should be downloaded. Supported values are:
#'
#'   \describe{
#'     \item{`both_sexes`}{Summary statistics generated using both male and
#'     female participants.}
#'     \item{`female`}{Summary statistics generated using female participants
#'     only.}
#'     \item{`male`}{Summary statistics generated using male participants
#'     only.}
#'   }
#'
#'   The default is:
#'
#'   \preformatted{
#'   c("both_sexes", "female", "male")
#'   }
#'
#'   which downloads summary statistics for all three sex groups. To download
#'   only one group, provide a single value, for example,
#'   `sex_groups = "both_sexes"`.
#'
#' @param dest_dir Character string. Directory where downloaded GWAS summary
#'   statistics files will be saved. The default is `"./"`, indicating that
#'   files will be saved in the current working directory.
#'
#' @return
#' Invisibly returns information about the downloaded files, if available.
#' The downloaded Neale Lab UK Biobank GWAS summary statistics files are
#' saved in the directory specified by `dest_dir`.
#'
#' @details
#' This function provides a convenient interface for downloading Neale Lab
#' UK Biobank GWAS summary statistics files for downstream polygenic risk
#' score (PRS) analysis.
#'
#' The phenotype code supplied through `pheno_code` corresponds to the
#' phenotype identifiers used in the UK Biobank GWAS Imputed v3 resource.
#' Users can refer to the UKBB GWAS Imputed v3 - File Manifest Release
#' 20180731 for:
#'
#' \itemize{
#'   \item phenotype codes;
#'   \item phenotype descriptions;
#'   \item available GWAS summary statistics files;
#'   \item annotations describing columns in the GWAS summary statistics
#'   files.
#' }
#'
#' The file manifest is available at:
#'
#' \url{https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU}
#'
#' The downloaded files contain GWAS summary statistics in bgzip-compressed
#' format (for example, `.bgz` or `.bgz.gz`) and can be used as input for
#' downstream PRS analyses, including [build_prs_weights()] and
#' [test_prs_association()].
#'
#' ## Examples of sex-specific downloads
#'
#' The default behavior:
#'
#' \preformatted{
#' download_neale_gwas(
#'   pheno_code = "100011_raw"
#' )
#' }
#'
#' downloads summary statistics for:
#'
#' \itemize{
#'   \item both sexes;
#'   \item females only;
#'   \item males only.
#' }
#'
#' To download only the GWAS summary statistics based on both sexes:
#'
#' \preformatted{
#' download_neale_gwas(
#'   pheno_code = "100011_raw",
#'   sex_groups = "both_sexes"
#' )
#' }
#'
#' To download female- and male-specific GWAS summary statistics:
#'
#' \preformatted{
#' download_neale_gwas(
#'   pheno_code = "100011_raw",
#'   sex_groups = c("female", "male")
#' )
#' }
#'
#' @examples
#' \dontrun{
#' ## Download GWAS summary statistics for all sex groups
#' download_neale_gwas(
#'   pheno_code = "100011_raw"
#' )
#'
#' ## Download both-sexes summary statistics only
#' download_neale_gwas(
#'   pheno_code = "100011_raw",
#'   sex_groups = "both_sexes"
#' )
#'
#' ## Save downloaded files to another directory
#' download_neale_gwas(
#'   pheno_code = "100011_raw",
#'   sex_groups = "both_sexes",
#'   dest_dir = "gwas_summary_statistics"
#' )
#' }
#'
#' @seealso
#' UKBB GWAS Imputed v3 - File Manifest Release 20180731:
#' \url{https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU},
#' [build_prs_weights()],
#' [test_prs_association()]
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




