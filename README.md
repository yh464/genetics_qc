## Genetics QC pipeline for UK Biobank
1. Use `ukb_extract.py` to extract a QC'ed list of subjects with the following quality control threshold: genetic sex agrees with reported sex; no excessive heterozygosity; European ancestry (self-reported and filtered for the first two genetic PCs within 5 SDs)
2. Use `genqc_snp.py` to filter SNPs using following criteria: MAF > 0.001 across UK Biobank, Hardy-Weinberg equilibrium exact test p-value (HWE) > 0.000001 across UK Biobank, and imputation INFO > 0.4.
3. Use `bfile_extract.py` to compile PLINK binaries for QC'ed subjects and filter SNPs using following criteria: MAF > 0.001 and HWE > 0.000001 across imaging cohort.
4. Imputed genotypes are obtained from UK Biobank.