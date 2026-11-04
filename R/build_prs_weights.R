
#' Build PRS SNP Weights from Source Neale Lab UK Biobank GWAS Results
#'
#' Extracts significant SNPs from a source Neale Lab UK Biobank GWAS summary
#' statistics file, annotates the selected SNPs using a Neale Lab-compatible
#' SNP annotation file, and generates a SNP weight file for downstream
#' polygenic risk score (PRS) calculation in target Neale Lab UK Biobank GWAS
#' datasets.
#'
#' The function applies a user-specified p-value threshold to identify
#' significant SNPs from the source GWAS summary statistics, retrieves SNP
#' annotation information, and saves SNP effect sizes as PRS weights. The
#' resulting weight file can subsequently be used to calculate PRS in target
#' GWAS datasets and evaluate associations between PRS and phenotypes.
#'
#' @param s_bgzfile Character string. Path to the compressed source Neale Lab
#'   UK Biobank GWAS summary statistics file in bgzip format (e.g.,
#'   `.bgz` or `.bgz.gz`). The file should contain SNP identifiers,
#'   effect alleles, effect size estimates, p-values, and other standard
#'   GWAS summary statistics.
#'
#' @param p_threshold Numeric value. P-value threshold used to select
#'   significant SNPs from the source GWAS summary statistics.
#'   For example, `5e-8` selects genome-wide significant variants.
#'
#' @param annfile Character string. Path to the SNP annotation file for
#'   Neale Lab UK Biobank GWAS summary statistics. The annotation file
#'   should contain SNP information required for PRS construction, including
#'   chromosome, genomic position, SNP identifier, and allele information.
#'
#' @param w_outfile Character string. Path to the output file where the
#'   selected SNPs and their corresponding effect-size weights will be
#'   written. The output file serves as the SNP weight file for downstream
#'   PRS calculation.
#' 
#' @param pop Character string. The reference population for identifying independent SNPs. 
#' Options are "EUR", "SAS", "EAS", "AFR", "AMR". 'legacy' 
#'
#' @return
#' A data frame containing the selected SNPs, annotation information, and
#' corresponding effect-size weights. The SNP weight file specified by
#' `w_outfile` is also written to disk.
#'
#' @details
#' ## OpenGWAS authentication
#'
#' This function uses the **ieugwasr** package to query the OpenGWAS API.
#' Therefore, a valid OpenGWAS JSON Web Token (JWT) must be stored in the
#' environment variable `OPENGWAS_JWT`.
#'
#' To check whether a JWT has been configured, run:
#'
#' \preformatted{
#' ieugwasr::get_opengwas_jwt()
#' }
#'
#' If no token is returned or the existing token has expired:
#'
#' \enumerate{
#'   \item Log in to the OpenGWAS website (https://api.opengwas.io) and generate a new JWT token.
#'   \item Open your `.Renviron` file by running:
#'
#'   \preformatted{
#'   usethis::edit_r_environ()
#'   }
#'
#'   \item Add or update the following line:
#'
#'   \preformatted{
#'   OPENGWAS_JWT="your_new_JWT_token"
#'   }
#'
#'   \item Save the file and restart your R session.
#' }
#'
#' ## Workflow
#'
#' The function performs the following steps:
#'
#' \enumerate{
#'   \item Reads the source Neale Lab UK Biobank GWAS summary statistics.
#'   \item Filters SNPs according to the specified p-value threshold.
#'   \item Matches the selected SNPs with the SNP annotation file.
#'   \item Extracts SNP effect sizes to construct PRS weights.
#'   \item Writes the SNP weight file specified by `w_outfile`.
#' }
#'
#' The generated SNP weight file can be used by downstream PRS functions to
#' calculate polygenic risk scores in target Neale Lab UK Biobank GWAS
#' datasets and test associations between PRS and target phenotypes.
#'
#' @examples
#' \dontrun{
#' ## Verify that an OpenGWAS JWT has been configured
#' ieugwasr::get_opengwas_jwt()
#'
#' ## Build a PRS weight file
#' build_prs_weights(
#'   s_bgzfile   = "source_gwas.txt.bgz",
#'   p_threshold = 5e-8,
#'   annfile     = "neale_snp_annotation.txt",
#'   w_outfile   = "prs_weights.txt"
#' )
#' }
#'
#' @seealso
#' ieugwasr::get_opengwas_jwt(),
#' usethis::edit_r_environ()
#'
#' @import ieugwasr 
#' @import data.table
#' 
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




