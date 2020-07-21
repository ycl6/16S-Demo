#' ---
#' title: "16S rDNA amplicon sequencing analysis using R (Part 4: picrust2 & ALDEx2)"
#' description: |
#'   Demonstration of 16S rRNA-based microbiome analysis using dada2, phyloseq, LEfSe, picrust2 and other tools
#' author:
#'   - name: "I-Hsuan Lin"
#'     url: https://github.com/ycl6
#'     affiliation: University of Manchester
#'     affiliation_url: https://www.manchester.ac.uk/
#' date: '`r format(Sys.Date(), "%B %d, %Y")`'
#' output:
#'     rmarkdown::html_document:
#'         theme: united
#'         highlight: tango
#'         self_contained: true
#'         toc: true
#'         toc_depth: 3
#'         toc_float:
#'             collapsed: false
#'             smooth_scroll: true
#' ---
#' 

#' 
#' -----
#' 
#' **16S-rDNA-V3-V4 repository:** https://github.com/ycl6/16S-rDNA-V3-V4
#' 
#' **picrust2:** [GitHub](https://github.com/picrust/picrust2), [Documentation](https://github.com/picrust/picrust2/wiki), [Paper](https://doi.org/10.1101/672295)
#' 
#' **ALDEx2**: [GitHub](https://github.com/ggloor/ALDEx_bioc), [Documentation](https://www.bioconductor.org/packages/release/bioc/html/ALDEx2.html), [Paper](https://doi.org/10.1186/2049-2618-2-15)
#' 
#' **Demo Dataset:** [PRJEB27564](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB27564) from *[Gut Microbiota in Parkinson's Disease: Temporal Stability and Relations to Disease Progression.](https://pubmed.ncbi.nlm.nih.gov/31221587/) EBioMedicine. 2019;44:691-707*
#' 
#' **License:** GPL-3.0
#' 
#' -----
#' 
#' # Workflow (Continues from Part 3)
#' 
#' > **Note:** `picrust2` requires Python 3.5 or 3.6
#' 
#' # Activate `picrust2` environment
#' 
#' > **Note:** Use `conda deactivate` to deactivate from a currently active conda environment if required
#' 
#' ```
#' conda activate picrust2
#' ```
#' 
#' # Run `picrust2`
#' 
#' * Requires fasta file `expr.asv.fasta` and biom file `expr.biom` generated in [Part 2: phyloseq](2_phyloseq_tutorial.html)
#' * The `--output` argument to specify the output folder for final files
#' * The `--processes N` argument to specify the number of CPUs to run picrust2 in parallel
#' 
#' ```
#' picrust2_pipeline.py --study_fasta outfiles/expr.asv.fasta --input outfiles/expr.biom \
#' 	--output picrust2_out_stratified --processes 6 \
#' 	--stratified --remove_intermediate --verbose
#' ```
#' 
#' > Completed PICRUSt2 pipeline in 2389.99 seconds.
#' 
#' > **Note:** Use `conda deactivate` to deactivate the `picrust2` environment if required
#' 
#' # Locate `picrust2` mapfiles
#' 
#' Use the `locate` command to locate the PATH that keeps the mapfiles
#' 
#' ```
#' locate description_mapfiles
#' ```
#' 
#' Example output:
#' 
#' ```
#' /home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles
#' /home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/KEGG_modules_info.tsv.gz
#' /home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/KEGG_pathways_info.tsv.gz
#' /home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/cog_info.tsv.gz
#' /home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/ec_level4_info.tsv.gz
#' /home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/ko_info.tsv.gz
#' /home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz
#' /home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/pfam_info.tsv.gz
#' /home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/tigrfam_info.tsv.gz
#' ```
#' 
#' # Start `R`
#' 
#' ```
#' cd /ngs/16S-Demo/PRJEB27564
#' 
#' R
#' ```
#' 
#' # Load package and set path
#' 
## ----load-libraries, message = FALSE------------------------------------------
library("data.table")	# Also requires R.utils to read gz and bz2 files
library("phyloseq")
library("ALDEx2")
library("dplyr")

#' 

#' 
#' ## Set `picrust2` output folder
#' 
#' Full path is not necessary if `R` is executed in the folder one-level above `picrust2_out_stratified`
#' 
## ----set-picrust2-output-folder, eval = FALSE---------------------------------
## picrust2 = "picrust2_out_stratified"

#' 
## -----------------------------------------------------------------------------
list.files(picrust2, recursive = TRUE)

#' 
#' ## Build `picrust2` output file paths
#' 
## ----set-picrust2-output-path-------------------------------------------------
p2_EC = paste0(picrust2, "/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_KO = paste0(picrust2, "/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_PW = paste0(picrust2, "/pathways_out/path_abun_unstrat.tsv.gz")

#' 
#' ## Set `picrust2` mapfile folder
#' 
#' Use the "description_mapfiles" PATH you located with the `locate` command above
#' 
## ----set-mapfile-folder, eval = FALSE-----------------------------------------
## mapfile = "/home/user/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles"

#' 
## -----------------------------------------------------------------------------
list.files(mapfile, recursive = TRUE)

#' 
#' ## Build `picrust2` map file paths
#' 
## ----set-mapfile-path---------------------------------------------------------
mapfile_EC = paste0(mapfile, "/ec_level4_info.tsv.gz")
mapfile_KO = paste0(mapfile, "/ko_info.tsv.gz")
mapfile_PW = paste0(mapfile, "/metacyc_pathways_info.txt.gz")

#' 
#' # Load saved workspace from Part 3
#' 
## ----load-image---------------------------------------------------------------
load("3_lefse_tutorial.RData")

#' 
#' # Prepare input data
#' 
#' ## Load map files 
#' 
## ----prepare-map--------------------------------------------------------------
mapEC = as.data.frame(fread(mapfile_EC, header = FALSE))
colnames(mapEC) = c("function","description")
mapKO = as.data.frame(fread(mapfile_KO, header = FALSE, sep = "\t"))
colnames(mapKO) = c("function","description")
mapPW = as.data.frame(fread(mapfile_PW, header = FALSE))
colnames(mapPW) = c("pathway","description")

#' 
#' ## Load `picrust2` output files
#' 
## ----prepare-input------------------------------------------------------------
p2EC = as.data.frame(fread(p2_EC))
rownames(p2EC) = p2EC$"function"
p2EC = as.matrix(p2EC[,-1])
p2EC = round(p2EC)

p2KO = as.data.frame(fread(p2_KO))
rownames(p2KO) = p2KO$"function"
p2KO = as.matrix(p2KO[,-1])
p2KO = round(p2KO)

p2PW = as.data.frame(fread(p2_PW))
rownames(p2PW) = p2PW$"pathway"
p2PW = as.matrix(p2PW[,-1])
p2PW = round(p2PW)

#' 
#' ## Subset results
#' 
#' * Subset-1: Patient vs. control at baseline
#' * Subset-2: Patient followup vs. patient baseline
#' 
## ----subset-results-----------------------------------------------------------
# Subset-1
p2EC1 = p2EC[,sample_names(ps1a)]
p2KO1 = p2KO[,sample_names(ps1a)]
p2PW1 = p2PW[,sample_names(ps1a)]

# Subset-2
p2EC2 = p2EC[,sample_names(ps1b)]
p2KO2 = p2KO[,sample_names(ps1b)]
p2PW2 = p2PW[,sample_names(ps1b)]

#' 
#' # Perform statistical analysis
#' 
#' We use the ANOVA-like differential expression (ALDEx2) compositional data analysis (CoDA) method to perform differential abundance testing **between 2 groups/conditions**
#' 
#' * If the test is ‘glm’, then effect should be FALSE. The ‘glm’ option evaluates the data as a one-way ANOVA using the glm and Kruskal-Wallace test
#' * If the test is ‘t’, then effect should be set to TRUE. The ‘t’ option evaluates the data as a two-factor experiment using both the Welch’s t and the Wilcoxon rank tests
#' * All tests include a Benjamini-Hochberg correction of the raw P values
#' 
#' The `mc.samples` argument defines the number of Monte Carlo samples to use to estimate the underlying distributions. The default is 128
#' 
#' **On Subset-1:** Patient vs. control at baseline
#' 
## ----run-ALDEx2-1, cache = TRUE-----------------------------------------------
set.seed(12345)
system.time({
        aldex2_EC1 = aldex(p2EC1, sample_data(ps1a)$Group, mc.samples = 500, test = "t", 
			   effect = TRUE, denom = "iqlr", verbose = TRUE)
})

set.seed(12345)
system.time({
        aldex2_KO1 = aldex(p2KO1, sample_data(ps1a)$Group, mc.samples = 500, test = "t", 
			   effect = TRUE, denom = "iqlr", verbose = TRUE)
})

set.seed(12345)
system.time({
        aldex2_PW1 = aldex(p2PW1, sample_data(ps1a)$Group, mc.samples = 500, test = "t", 
			   effect = TRUE, denom = "iqlr", verbose = TRUE)
})

#' 
#' **On Subset-2:** Patient followup vs. patient baseline
#' 
## ----run-ALDEx2-2, cache = TRUE-----------------------------------------------
set.seed(12345)
system.time({
        aldex2_EC2 = aldex(p2EC2, sample_data(ps1b)$Time, mc.samples = 500, test = "t", 
			   effect = TRUE, denom = "iqlr", verbose = TRUE)
})

set.seed(12345)
system.time({
        aldex2_KO2 = aldex(p2KO2, sample_data(ps1b)$Time, mc.samples = 500, test = "t", 
			   effect = TRUE, denom = "iqlr", verbose = TRUE)
})

set.seed(12345)
system.time({
        aldex2_PW2 = aldex(p2PW2, sample_data(ps1b)$Time, mc.samples = 500, test = "t", 
			   effect = TRUE, denom = "iqlr", verbose = TRUE)
})

#' 
#' ## See a ALDEx2 output
#' 
## -----------------------------------------------------------------------------
head(aldex2_EC1)

#' 
#' ## Check estimated effect size
#' 
#' > ALDEx2 authors suggest that an effect size of 1 or greater can be used as significance cutoff
#' 
#' **On Subset-1:** Patient vs. control at baseline
#' 
#' > None of the metabolic predictions show differential abundances between Parkinson's patients and control subjects at baseline
#' 
## -----------------------------------------------------------------------------
quantile(aldex2_EC1$effect, seq(0, 1, 0.1))

#' 
## -----------------------------------------------------------------------------
quantile(aldex2_KO1$effect, seq(0, 1, 0.1))

#' 
## -----------------------------------------------------------------------------
quantile(aldex2_PW1$effect, seq(0, 1, 0.1))

#' 
#' **On Subset-2:** Patient followup vs. patient baseline
#' 
#' > None of the metabolic predictions show differential abundances between Parkinson's patients at baseline and at followup
#' 
## -----------------------------------------------------------------------------
quantile(aldex2_EC2$effect, seq(0, 1, 0.1))

#' 
## -----------------------------------------------------------------------------
quantile(aldex2_KO2$effect, seq(0, 1, 0.1))

#' 
## -----------------------------------------------------------------------------
quantile(aldex2_PW2$effect, seq(0, 1, 0.1))

#' 
#' # Plotting of outputs
#' 
#' ## Create MW and MA plots
#' 
#' Use `aldex.plot` function to create MW (fold-change to variance/effect) and MA (Bland-Altman) plots
#' 
#' The plot below shows the relationship between effect size and P values and BH-adjusted P values in the test dataset.
#' 
#' **On Subset-1:** Patient vs. control at baseline
#' 
## ----plot-MW-MA-1-------------------------------------------------------------
png("images/ALDEx2_picrust2_MW_MA_1.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow=c(3,2))
aldex.plot(aldex2_EC1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
	   called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(EC) MW Plot")

aldex.plot(aldex2_EC1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
	   called.cex = 0.6, cutoff = 0.05, xlab = "Log-ratio abundance", ylab = "Difference")
title(main = "(EC) MA Plot")

aldex.plot(aldex2_KO1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
	   called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(KO) MW Plot")

aldex.plot(aldex2_KO1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
	   called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(KO) MA Plot")

aldex.plot(aldex2_PW1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
	   called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(PW) MW Plot")

aldex.plot(aldex2_PW1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
	   called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(PW) MA Plot")
invisible(dev.off())

#' 

#' 
#' **On Subset-2:** Patient followup vs. patient baseline
#' 
## ----plot-MW-MA-2-------------------------------------------------------------
png("images/ALDEx2_picrust2_MW_MA_2.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow=c(3,2))
aldex.plot(aldex2_EC2, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
	   called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(EC) MW Plot")

aldex.plot(aldex2_EC2, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
	   called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(EC) MA Plot")

aldex.plot(aldex2_KO2, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
	   called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(KO) MW Plot")

aldex.plot(aldex2_KO2, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
	   called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(KO) MA Plot")

aldex.plot(aldex2_PW2, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
	   called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(PW) MW Plot")

aldex.plot(aldex2_PW2, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
	   called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(PW) MA Plot")
invisible(dev.off())

#' 

#' 
#' ## Relationship between effect, difference, and P values
#' 
#' **On Subset-1:** Patient vs. control at baseline
#' 
## ----plot-P-adjP-1------------------------------------------------------------
png("images/ALDEx2_picrust2_P_adjP_1.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow=c(3,2))
plot(aldex2_EC1$effect, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(EC) Effect size plot")
points(aldex2_EC1$effect, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_EC1$diff.btw, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(EC) Volcano plot")
points(aldex2_EC1$diff.btw, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_KO1$effect, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(KO) Effect size plot")
points(aldex2_KO1$effect, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_KO1$diff.btw, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(KO) Volcano plot")
points(aldex2_KO1$diff.btw, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_PW1$effect, aldex2_PW1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(PW) Effect size plot")
points(aldex2_PW1$effect, aldex2_PW1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_PW1$diff.btw, aldex2_PW1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(PW) Volcano plot")
points(aldex2_PW1$diff.btw, aldex2_PW1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
invisible(dev.off())

#' 

#' 
#' **On Subset-2:** Patient followup vs. patient baseline
#' 
## ----plot-P-adjP-2------------------------------------------------------------
png("images/ALDEx2_picrust2_P_adjP_2.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow=c(3,2))
plot(aldex2_EC2$effect, aldex2_EC2$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(EC) Effect size plot")
points(aldex2_EC2$effect, aldex2_EC2$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_EC2$diff.btw, aldex2_EC2$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(EC) Volcano plot")
points(aldex2_EC2$diff.btw, aldex2_EC2$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_KO2$effect, aldex2_KO2$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(KO) Effect size plot")
points(aldex2_KO2$effect, aldex2_KO2$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_KO2$diff.btw, aldex2_KO2$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(KO) Volcano plot")
points(aldex2_KO2$diff.btw, aldex2_KO2$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_PW2$effect, aldex2_PW2$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(PW) Effect size plot")
points(aldex2_PW2$effect, aldex2_PW2$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_PW2$diff.btw, aldex2_PW2$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(PW) Volcano plot")
points(aldex2_PW2$diff.btw, aldex2_PW2$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
invisible(dev.off())

#' 

#' 
#' # Create output files
#' 
#' ## Merge with map file data
#' 
## ----merge-files--------------------------------------------------------------
# Subset-1
df_EC1 = aldex2_EC1 %>% tibble::rownames_to_column(var = "EC") %>% 
	inner_join(mapEC, by = c("EC" = "function")) %>% arrange(EC)

df_KO1 = aldex2_KO1 %>% tibble::rownames_to_column(var = "KO") %>% 
	inner_join(mapKO, by = c("KO" = "function")) %>% arrange(KO)

df_PW1 = aldex2_PW1 %>% tibble::rownames_to_column(var = "Pathway") %>% 
	inner_join(mapPW, by = c("Pathway" = "pathway")) %>% arrange(Pathway)

# Subset-2
df_EC2 = aldex2_EC2 %>% tibble::rownames_to_column(var = "EC") %>%
        inner_join(mapEC, by = c("EC" = "function")) %>% arrange(EC)

df_KO2 = aldex2_KO2 %>% tibble::rownames_to_column(var = "KO") %>%
        inner_join(mapKO, by = c("KO" = "function")) %>% arrange(KO)

df_PW2 = aldex2_PW2 %>% tibble::rownames_to_column(var = "Pathway") %>%
        inner_join(mapPW, by = c("Pathway" = "pathway")) %>% arrange(Pathway)

#' 
#' ## Output to file
#' 
## ----output-files-------------------------------------------------------------
# Subset-1
write.table(df_EC1, file = "outfiles/ALDEx2_picrust2_EC_results_1.txt", sep = "\t", quote = F, 
	    row.names = F, col.names = T)
write.table(df_KO1, file = "outfiles/ALDEx2_picrust2_KO_results_1.txt", sep = "\t", quote = F, 
	    row.names = F, col.names = T)
write.table(df_PW1, file = "outfiles/ALDEx2_picrust2_Pathway_results_1.txt", sep = "\t", quote = F, 
	    row.names = F, col.names = T)

# Subset-2
write.table(df_EC2, file = "outfiles/ALDEx2_picrust2_EC_results_2.txt", sep = "\t", quote = F, 
	    row.names = F, col.names = T)
write.table(df_KO2, file = "outfiles/ALDEx2_picrust2_KO_results_2.txt", sep = "\t", quote = F, 
	    row.names = F, col.names = T)
write.table(df_PW2, file = "outfiles/ALDEx2_picrust2_Pathway_results_2.txt", sep = "\t", quote = F, 
	    row.names = F, col.names = T)

#' 
#' # Save current workspace
#' 
## ----save-image---------------------------------------------------------------
save.image(file = "4_picrust2_tutorial.RData")

#' 
#' # Session information
#' 
## ----session-info-------------------------------------------------------------
sessionInfo()

