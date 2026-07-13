# Test script for cgiarGenomics::impute_beagle function

# Paths
vcf_path <- "tests/vcf_fmt/diploid.vcf.gz"
file_imputed_path <- "file40bc782a727c_imputed.vcf.gz"
jre_path <- "C:/Program Files/Java/jre-1.8"
beagle_path <- "./beagle.27Feb25.75f.jar"

cat("Reading input VCF with read_vcf...\n")
gl_input <- read_vcf(vcf_path, ploidity = 2, na_reps = "./.")
print(gl_input)

cat("\nImputing with cgiarGenomics::impute_beagle...\n")
gl_imputed_func <- impute_beagle(gl_input,
                                 jre_path = jre_path,
                                 beagle_path = beagle_path,
                                 memory = "Xmx5g",
                                 nthreads = 4) # use 16 threads to match manual run

gl <- read_vcf("~/beagle_out_40bc233a468c.vcf/beagle_out_40bc233a468c.vcf", sep = "\\|")
print(gl_imputed_func)

cat("\nReading the manual Beagle imputed VCF file...\n")
# Replace | with / on a copy of the manual imputed VCF so read_vcf can parse it
temp_orig_imputed <- tempfile(fileext = ".vcf.gz")
con <- gzfile(file_imputed_path, "rt")
lines <- readLines(con)
close(con)
lines_replaced <- gsub("|", "/", lines, fixed = TRUE)
con_out <- gzfile(temp_orig_imputed, "wt")
writeLines(lines_replaced, con_out)
close(con_out)

gl_imputed_file <- read_vcf(temp_orig_imputed, ploidity = 2, na_reps = "./.")
print(gl_imputed_file)
unlink(temp_orig_imputed)

# Comparison
cat("\n--- Comparison Results ---\n")
m_func <- as.matrix(gl_imputed_func)
m_file <- as.matrix(gl_imputed_file)

match_matrix <- all(m_func == m_file, na.rm = TRUE) && identical(is.na(m_func), is.na(m_file))
cat("1. Do genotype matrices match exactly?", match_matrix, "\n")

match_chrom <- identical(adegenet::chromosome(gl_imputed_func), adegenet::chromosome(gl_imputed_file))
cat("2. Do chromosomes match exactly?", match_chrom, "\n")

match_pos <- identical(adegenet::position(gl_imputed_func), adegenet::position(gl_imputed_file))
cat("3. Do positions match exactly?", match_pos, "\n")

match_alleles <- identical(adegenet::alleles(gl_imputed_func), adegenet::alleles(gl_imputed_file))
cat("4. Do alleles match exactly?", match_alleles, "\n")

if (match_matrix && match_chrom && match_pos && match_alleles) {
  cat("\nSUCCESS: The output of impute_beagle matches file40bc782a727c_imputed.vcf.gz perfectly!\n")
} else {
  cat("\nFAILURE: Mismatch found between function output and manual imputation file.\n")
}
