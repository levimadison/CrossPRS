
#' Build PRS SNP Weights from Neale Lab UK Biobank GWAS Results
#'
#' Extracts significant SNPs from a source Neale Lab UK Biobank GWAS summary
#' statistics file, annotates the selected SNPs using a Neale Lab variant
#' annotation file, and writes a SNP weight file for downstream polygenic
#' risk score (PRS) analysis.
#'
#' The source GWAS file can be obtained using [download_neale_gwas()] or
#' downloaded directly from the Neale Lab UK Biobank GWAS Imputed v3
#' resources. The annotation file is a separate input.
#'
#' @param s_bgzfile Character string. Path to a downloaded, bgzip-compressed
#'   source Neale Lab UK Biobank GWAS summary statistics file, such as a
#'   \code{.tsv.bgz} file. The file should contain the summary-statistic
#'   columns required by this function. Use [download_neale_gwas()] to
#'   download a phenotype file, or find its download link in the
#'   \code{Manifest 201807} sheet described in Details. Supply the downloaded
#'   file's local path to this argument.
#'
#' @param p_threshold Numeric value. P-value threshold used to select
#'   significant SNPs from the source GWAS summary statistics. For example,
#'   \code{5e-8} selects variants with p-values below the genome-wide
#'   significance threshold.
#'
#' @param annfile Character string. Path to the downloaded Neale Lab variant
#'   annotation file. The file should contain the variant and allele
#'   information required by this function. The Neale Lab annotation file is
#'   named \code{variants.tsv.bgz} and is available at:
#'   \url{https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz}
#'
#'   Download the annotation file separately and supply its local path to
#'   this argument. This file is not the GWAS summary-statistics file
#'   downloaded by [download_neale_gwas()].
#'
#' @param w_outfile Character string. Path to the output file (text file) where the
#'   selected SNPs, annotation information, and corresponding GWAS
#'   summary statistics will be written. This file can be used by downstream PRS
#'   functions.
#'
#' @param pop Character string. Reference population used to identify
#'   independent SNPs. Supported values are \code{"EUR"}, \code{"SAS"},
#'   \code{"EAS"}, \code{"AFR"},  \code{"AMR"} and \code{"legacy"}. 
#'
#' @return A data frame containing the selected SNPs, annotation information,
#'   and corresponding GWAS summary statistics. The SNP weight file specified by
#'   \code{w_outfile} is also written to disk.
#'
#' @details
#' ## Neale Lab input files
#'
#' The Neale Lab UK Biobank GWAS Imputed v3 File Manifest Release 20180731
#' provides information about available GWAS files. Its \code{README} sheet
#' describes the files, and its \code{Manifest 201807} sheet lists phenotype
#' codes, phenotype descriptions, and download links:
#'
#' \url{https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU}
#'
#' For example, female-specific body mass index (BMI) summary-statistics
#' files include:
#' \itemize{
#'   \item \code{21001_irnt}: inverse rank normalized BMI:
#'     \url{https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_irnt.gwas.imputed_v3.female.tsv.bgz}
#'   \item \code{21001_raw}: untransformed BMI:
#'     \url{https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_raw.gwas.imputed_v3.female.tsv.bgz}
#' }
#'
#' These are examples of direct download links. Use the manifest to find the
#' file for the phenotype and sex group required by your analysis. Alternatively,
#' use [download_neale_gwas()] to download the GWAS summary-statistics file,
#' then pass its local path as \code{s_bgzfile}.
#'
#' The variant annotation file is downloaded separately. The Neale Lab
#' \code{variants.tsv.bgz} file is available at:
#' \url{https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz}
#'
#' ## OpenGWAS authentication
#'
#' This function uses the \pkg{ieugwasr} package to query the OpenGWAS API.
#' If the function's OpenGWAS query requires authentication, a valid OpenGWAS
#' JSON Web Token (JWT) must be available in the \code{OPENGWAS_JWT}
#' environment variable.
#'
#' To check whether a JWT is configured, run:
#'
#' \preformatted{
#' ieugwasr::get_opengwas_jwt()
#' }
#'
#' If a token is missing or expired:
#' \enumerate{
#'   \item Obtain a token from \url{https://api.opengwas.io}.
#'   \item Open your \code{.Renviron} file by running:
#'     \preformatted{
#'     usethis::edit_r_environ()
#'     }
#'   \item Add or update the following line:
#'     \preformatted{
#'     OPENGWAS_JWT="your_new_JWT_token"
#'     }
#'   \item Save the file and restart your R session.
#' }
#'
#' Do not share or commit your JWT.
#'
#' ## Workflow
#'
#' The function:
#' \enumerate{
#'   \item Reads the source GWAS summary statistics from \code{s_bgzfile}.
#'   \item Selects SNPs according to \code{p_threshold}.
#'   \item Matches the selected SNPs to the annotation file in \code{annfile}.
#'   \item Extracts the information needed to construct PRS SNP weights.
#'   \item Writes the resulting weight file to \code{w_outfile}.
#' }
#'
#' The output weight file can be used by downstream PRS functions, such as
#' [test_prs_association()], to test associations in a target GWAS.
#'
#' @examples
#' \dontrun{
#' ## Download a source GWAS file with download_neale_gwas(), or download it
#' ## directly using a link from the Neale Lab manifest.
#' download_neale_gwas(
#'   pheno_code = "21001_irnt",
#'   sex_groups = "female",
#'   dest_dir = "gwas_data"
#' )
#'
#' ## Download variants.tsv.bgz separately from the Neale Lab annotation URL
#' ## documented above, then provide both local file paths.
#' build_prs_weights(
#'   s_bgzfile = "gwas_data/21001_irnt.gwas.imputed_v3.female.tsv.bgz",
#'   p_threshold = 5e-8,
#'   annfile = "gwas_data/variants.tsv.bgz",
#'   w_outfile = "prs_weights.txt",
#'   pop = "EUR"
#' )
#' }
#'
#' @seealso
#' [download_neale_gwas()],
#' [test_prs_association()],
#' \code{ieugwasr::get_opengwas_jwt()},
#' \code{usethis::edit_r_environ()}
#'
#' @import ieugwasr
#' @import data.table
#' @export
build_prs_weights<-function(s_bgzfile, p_threshold = 5e-8,annfile="variants.tsv.bgz", w_outfile, pop = "EUR") {
	#library(data.table)
	#library(ieugwasr)
	
# -------------------------------------------------------------------------
# Read and Filter GWAS Summary Stats
# -------------------------------------------------------------------------
	
	cat("Reading and filtering GWAS summary statistics...\n")
	
	gwas_data <- data.table::fread(
			s_bgzfile,
			select = c("variant", "minor_AF", "beta", "se", "pval"), 
			colClasses = c(variant = "character", pval = "numeric")
	)
	
# Keep variants with p-value <= p_threshold
	filtered_gwas <- gwas_data[pval <= p_threshold]
	
	
# -------------------------------------------------------------------------
# Parse "chr:pos:ref:alt" into Dedicated Columns
# -------------------------------------------------------------------------
	cat("Parsing variant coordinates...\n")
	filtered_gwas[, c("chr", "pos", "ref", "alt") := tstrsplit(variant, ":", fixed = TRUE)]
	filtered_gwas[, pos := as.integer(pos)]
	
	
# -------------------------------------------------------------------------
# Map variant (e.g. "1:15791:C:T") to "rsid"
# -------------------------------------------------------------------------
	cat("Loading variant manifest to map coordinates to rsIDs...\n")
	manifest <- data.table::fread(annfile, select = c("variant", "rsid"))
	
# Merge the rsIDs using the coordinate string as the key
	filtered_gwas_mapped <- merge(filtered_gwas, manifest, by = "variant", all.x = TRUE)
	
# Remove any variants that lack a valid rsID mapping (missing from dbSNP)
	filtered_gwas_mapped <- filtered_gwas_mapped[!is.na(rsid) & rsid != ""]
	
# -------------------------------------------------------------------------
# Run ld_clump with True rsIDs
# -------------------------------------------------------------------------
	cat("Preparing clean input for OpenGWAS API...\n")
	
# The 'rsid' column contains real rsIDs (e.g., rs12345),
# and the 'pval' column contains the p-values.
	clump_input <- data.frame(
			rsid = filtered_gwas_mapped$rsid, 
			pval = filtered_gwas_mapped$pval,
			stringsAsFactors = FALSE
	)
	
	cat("Performing remote LD clumping...\n")
	independent_snps <- ieugwasr::ld_clump(
			dat = clump_input,
			clump_kb = 10000,   
			clump_r2 = 0.001,   
			clump_p = 1,        
			pop = pop        
	)
	
# -------------------------------------------------------------------------
# Isolate  final independent variants
# -------------------------------------------------------------------------
# Map back to get final list of independent variants with all details
	final_prs_variants <- filtered_gwas_mapped[rsid %in% independent_snps$rsid]
	
	print(head(final_prs_variants))
	fwrite(final_prs_variants, w_outfile, sep = "\t")
}




