#' ---
#' title: "16S rDNA amplicon sequencing analysis using R (Part 2: phyloseq)"
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
#' **phyloseq:** [GitHub](https://github.com/joey711/phyloseq), [Documentation](https://joey711.github.io/phyloseq/index.html), [Paper](https://doi.org/10.1371/journal.pone.0061217)
#' 
#' **GUniFrac:** [Paper](https://doi.org/10.1093/bioinformatics/bts342)
#' 
#' **Demo Dataset:** [PRJEB27564](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB27564) from *[Gut Microbiota in Parkinson's Disease: Temporal Stability and Relations to Disease Progression.](https://pubmed.ncbi.nlm.nih.gov/31221587/) EBioMedicine. 2019;44:691-707*
#' 
#' **License:** GPL-3.0
#' 
#' -----
#' 
#' # Workflow (Continues from Part 1)
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
library("data.table")
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("ggbeeswarm")
library("ggrepel")
library("vegan")
library("tidyverse")

#' 
#' Set the full-path to additional scripts
#' 
#' > **Note:** Remember to change the PATH to the script accordingly
#' 
## ----set-taxa_summary-Rscript-path, eval = FALSE------------------------------
## source("/path-to-script/taxa_summary.R", local = TRUE)

#' 

#' 
#' # Load saved workspace from Part 1
#' 
## ----load-image---------------------------------------------------------------
load("1_dada2_tutorial.RData")

#' 
#' Check `phyloseq` object
#' 
## -----------------------------------------------------------------------------
ps

#' 
#' # Filter Taxa
#' 
#' ## Subset `ps` using `subset_taxa`
#' 
#' Here, we limit the taxa to those from Bacteria, Phylum and Class in not unassigned, and not belonging to the Mitochondria Family
#' 
## ----subset_taxa--------------------------------------------------------------
ps0 = subset_taxa(ps, Kingdom == "Bacteria" & !is.na(Phylum) & !is.na(Class) & 
		  Family != "Mitochondria")

#' 
#' ## Create a total counts data.table
#' 
## ----create-tdt, warning = FALSE, fig.width = 6, fig.height = 6, fig.align = "center", dpi = 100----
tdt = data.table(setDT(as.data.frame(tax_table(ps0))), 
		 TotalCounts = taxa_sums(ps0), SV = taxa_names(ps0))

tdt

ggplot(tdt, aes(TotalCounts)) + geom_histogram(bins = 50) + theme_bw() + 
	ggtitle("Histogram of Total Counts")

tdt[(TotalCounts == 1), .N]     # singletons

tdt[(TotalCounts == 2), .N]     # doubletons

# taxa cumulative sum
taxcumsum = tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]

# Plot the cumulative sum of ASVs against the total counts
pCumSum = ggplot(taxcumsum, aes(TotalCounts, CumSum)) + geom_point() + theme_bw() + 
	xlab("Filtering Threshold") + ylab("ASV Filtered")

png("images/FilterTaxa-taxa-abundance.png", width = 8, height = 8, units = "in", res = 300)
gridExtra::grid.arrange(pCumSum, pCumSum + xlim(0, 500), 
			pCumSum + xlim(0, 100), pCumSum + xlim(0, 50), nrow = 2, 
			top = "ASVs that would be filtered vs. minimum taxa counts threshold")
invisible(dev.off())

#' 

#' 
#' ## Create a prevalence data.table
#' 
## ----create-mdt, warning = FALSE, fig.width = 6, fig.height = 6, fig.align = "center", dpi = 100----
mdt = fast_melt(ps0)

mdt

prevdt = mdt[, list(Prevalence = sum(count > 0), TotalCounts = sum(count), 
		    MaxCounts = max(count)), by = TaxaID]

prevdt

ggplot(prevdt, aes(Prevalence)) + geom_histogram(bins = 50) + theme_bw() +
	ggtitle("Histogram of Taxa Prevalence")

prevdt[(Prevalence == 1), .N]	# singletons

prevdt[(Prevalence == 2), .N]	# doubletons

ggplot(prevdt, aes(MaxCounts)) + geom_histogram(bins = 100) + xlim(0, 500) + theme_bw() + 
	ggtitle("Histogram of Maximum TotalCounts")

table(prevdt$MaxCounts)[1:50]

prevdt[(MaxCounts == 1)]	# singletons

# taxa cumulative sum
prevcumsum = prevdt[, .N, by = Prevalence]
setkey(prevcumsum, Prevalence)
prevcumsum[, CumSum := cumsum(N)]

# Plot the cumulative sum of ASVs against the prevalence
pPrevCumSum = ggplot(prevcumsum, aes(Prevalence, CumSum)) + geom_point(size = 2, alpha = 0.5) + 
	theme_bw() + xlab("Filtering Threshold") + ylab("ASVs Filtered") + 
	ggtitle("ASVs that would be filtered vs. minimum sample count threshold")

png("images/FilterTaxa-taxa-prevalence.png", width = 8, height = 8, units = "in", res = 300)
pPrevCumSum
invisible(dev.off())

#' 

#' 
#' Prevalence vs. Total Count Scatter plot
#' 
## ----warning = FALSE----------------------------------------------------------
png("images/FilterTaxa-Prevalence-TotalCounts.png", width = 8, height = 8, units = "in", res = 300)
ggplot(prevdt, aes(Prevalence, TotalCounts)) + geom_point(size = 2, alpha = 0.5) + 
	scale_y_log10() + theme_bw() + xlab("Prevalence [No. Samples]") + ylab("TotalCounts [Taxa]")
invisible(dev.off())

#' 

#' 
#' Colored by phylum
#' 
## ----addPhylum, warning = FALSE-----------------------------------------------
addPhylum = unique(copy(mdt[, list(TaxaID, Phylum)]))

# Join by TaxaID
setkey(prevdt, TaxaID)
setkey(addPhylum, TaxaID)
prevdt = addPhylum[prevdt]
setkey(prevdt, Phylum)

png("images/FilterTaxa-Prevalence-TotalCounts-Phylum.png", width = 8, height = 8, 
    units = "in", res = 300)
ggplot(prevdt, aes(Prevalence, TotalCounts, color = Phylum)) + 
	geom_point(size = 1, alpha = 0.5) + scale_y_log10() + theme_bw() + 
	facet_wrap(~Phylum, nrow = 4) + theme(legend.position="none") + 
	xlab("Prevalence [No. Samples]") + ylab("Total Abundance")
invisible(dev.off())

#' 

#' 
#' ## Define filter threshold
#' 
#' > Remember to review the **images/FilterTaxa-\*.png** plots to select the best paramters
#' 
## ----filter-threshold---------------------------------------------------------
prevalenceThreshold = 5
abundanceThreshold = 10
maxThreshold = 5

#' 
#' ## Perform taxa filtering
#' 
## ----filter-ps0---------------------------------------------------------------
keepTaxa = prevdt[(Prevalence > prevalenceThreshold & 
		   TotalCounts > abundanceThreshold & 
		   MaxCounts > maxThreshold), TaxaID]

ps1 = prune_taxa(keepTaxa, ps0)

phy_tree(ps1) = ape::root(phy_tree(ps1), sample(taxa_names(ps1), 1), resolve.root = TRUE)

ps1

#' 
#' # Create output files
#' 
#' ## Create FASTA file
#' 
## ----output-fasta, warning = FALSE, cache = TRUE------------------------------
dada2::uniquesToFasta(df[rownames(df) %in% taxa_names(ps1),], "outfiles/expr.asv.fasta", 
		      ids = df[rownames(df) %in% taxa_names(ps1),]$id)

#' 
#' ## Create BIOM file
#' 
## ----output-biom, warning = FALSE, cache = TRUE-------------------------------
biomformat::write_biom(biomformat::make_biom(data = t(as.matrix(otu_table(ps1)))), 
		       "outfiles/expr.biom")

#' 
#' ## Create ASV and Taxonomy tables
#' 
## ----output-otu-tax-tables, warning = FALSE, cache = TRUE---------------------
write.table(as.data.table(otu_table(ps1), keep.rownames = T), file = "outfiles/expr.otu_table.txt", 
	    sep = "\t", quote = F, row.names = F, col.names = T)
write.table(as.data.table(tax_table(ps1), keep.rownames = T), file = "outfiles/expr.tax_table.txt", 
	    sep = "\t", quote = F, row.names = F, col.names = T)

#' 
#' ## Create raw abundance tables
#' 
#' 
## ----output-abund-tables, warning = FALSE, cache = TRUE-----------------------
write.table(ps1 %>% psmelt() %>% arrange(OTU) %>% rename(ASV = OTU) %>%
	    select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species, Sample_ID, Abundance) %>%
	    spread(Sample_ID, Abundance),
	    file = "outfiles/expr.abundance.all.txt", sep = "\t", 
	    quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Phylum") %>% psmelt() %>%
	    select(Phylum, Sample_ID, Abundance) %>% spread(Sample_ID, Abundance),
	    file = "outfiles/expr.abundance.abphy.txt", sep = "\t", 
	    quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Class") %>% psmelt() %>%
	    select(Class, Sample_ID, Abundance) %>% spread(Sample_ID, Abundance),
	    file = "outfiles/expr.abundance.abcls.txt", sep = "\t", 
	    quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Family") %>% psmelt() %>%
	    select(Family, Sample_ID, Abundance) %>% spread(Sample_ID, Abundance),
	    file = "outfiles/expr.abundance.abfam.txt", sep = "\t", 
	    quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Genus") %>% psmelt() %>%
	    select(Genus, Sample_ID, Abundance) %>% spread(Sample_ID, Abundance),
	    file = "outfiles/expr.abundance.abgen.txt", sep = "\t", 
	    quote = F, row.names = F, col.names = T)

#' 
#' ## Create relative abundance tables
#' 
#' Use the `transform_sample_counts` function to transforms the sample counts of a taxa abundance matrix according to the provided function, in the case the fractions of the whole sum
#' 
## ----output-rel-abund-tables, warning = FALSE, cache = TRUE-------------------
write.table(ps1 %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>%
	    arrange(OTU) %>% rename(ASV = OTU) %>%
	    select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species, Sample_ID, Abundance) %>%
	    spread(Sample_ID, Abundance),
    file = "outfiles/expr.relative_abundance.all.txt", 
    sep = "\t", quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Phylum") %>%
	    transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>%
	    select(Phylum, Sample_ID, Abundance) %>% spread(Sample_ID, Abundance),
    file = "outfiles/expr.relative_abundance.abphy.txt", 
    sep = "\t", quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Class") %>%
	    transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>%
	    select(Class, Sample_ID, Abundance) %>% spread(Sample_ID, Abundance),
    file = "outfiles/expr.relative_abundance.abcls.txt", 
    sep = "\t", quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Family") %>%
	    transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>%
	    select(Family, Sample_ID, Abundance) %>% spread(Sample_ID, Abundance),
    file = "outfiles/expr.relative_abundance.abfam.txt", 
    sep = "\t", quote = F, row.names = F, col.names = T)
write.table(ps1 %>% tax_glom(taxrank = "Genus") %>%
	    transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>%
	    select(Genus, Sample_ID, Abundance) %>% spread(Sample_ID, Abundance),
    file = "outfiles/expr.relative_abundance.abgen.txt", 
    sep = "\t", quote = F, row.names = F, col.names = T)

#' 
#' # Plot abundance
#' 
## ----prepare-plot-bar---------------------------------------------------------
# Transform to proportions (relative abundances)
ps1.rp = transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))

# Top N taxa
N = 200
topN = names(sort(taxa_sums(ps1), decreasing=TRUE))[1:N]
ps1.topN = transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps1.topN = prune_taxa(topN, ps1.topN)

ps1.topN

#' 
#' ## At Phylum level
#' 
## ----plot-bar-phylum, cache = TRUE--------------------------------------------
# All taxa
ptaxa1 = plot_bar(ps1.rp, x = "Sample_ID", fill = "Phylum", 
		  title = paste(ntaxa(ps1.rp), "Taxa (Phylum)")) + 
	geom_bar(stat = "identity", size = 0.1, color = "black") + 
	facet_wrap(~Group2, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) + 
	scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
	theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7), 
	      axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1), 
	      axis.text.y = element_text(size = 12)) + xlab("Sample") + ylab("Relative abundance")

# Top 200 taxa
ptaxa2 = plot_bar(ps1.topN, x = "Sample_ID", fill = "Phylum", 
		  title = paste("Top",N, "Taxa (Phylum)")) + 
	geom_bar(stat = "identity", size = 0.1, color = "black") + 
	facet_wrap(~Group2, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) + 
	scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
	theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7), 
	      axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1), 
	      axis.text.y = element_text(size = 12)) + xlab("Sample") + ylab("Relative Abundance")

png("images/plot_bar_phylum.png", width = 10, height = 8, units = "in", res = 300)
gridExtra::grid.arrange(ptaxa1, ptaxa2, nrow = 2)
invisible(dev.off())

#' 

#' 
#' ## At Class level
#' 
## ----plot-bar-class, cache = TRUE---------------------------------------------
ptaxa3 = plot_bar(ps1.rp, x = "Sample_ID", fill = "Class", 
		  title = paste(ntaxa(ps1.rp), "Taxa (Class)")) +
        geom_bar(stat = "identity", size = 0.1, color = "black") +
        facet_wrap(~Group2, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) +
        scale_y_continuous(breaks = seq(0, 1, 0.1)) +
        theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7),
              axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
              axis.text.y = element_text(size = 12)) + xlab("Sample") + ylab("Relative abundance")

ptaxa4 = plot_bar(ps1.topN, x = "Sample_ID", fill = "Class", 
		  title = paste("Top",N, "Taxa (Class)")) +
        geom_bar(stat = "identity", size = 0.1, color = "black") +
        facet_wrap(~Group2, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) +
        scale_y_continuous(breaks = seq(0, 1, 0.1)) +
        theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7),
              axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
              axis.text.y = element_text(size = 12)) + xlab("Sample") + ylab("Relative Abundance")

png("images/plot_bar_class.png", width = 10, height = 8, units = "in", res = 300)
gridExtra::grid.arrange(ptaxa3, ptaxa4, nrow = 2)
invisible(dev.off())

#' 

#' 
#' ## At Family level
#' 
## ----plot-bar-family, cache = TRUE--------------------------------------------
ptaxa5 = plot_bar(ps1.topN, x = "Sample_ID", fill = "Family", 
		  title = paste("Top",N, "Taxa (Family)")) +
        geom_bar(stat = "identity", size = 0.1, color = "black") +
        facet_wrap(~Group2, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 1)) +
        scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
        theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7),
              axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
              axis.text.y = element_text(size = 12)) + xlab("Sample") + ylab("Relative Abundance")

png("images/plot_bar_family.png", width = 10, height = 8, units = "in", res = 300)
ptaxa5
invisible(dev.off())

#' 

#' 
#' ## At Genus level
#' 
## ----plot-bar-genus, cache = TRUE---------------------------------------------
ptaxa6 = plot_bar(ps1.topN, x = "Sample_ID", fill = "Genus", 
		  title = paste("Top",N, "Taxa (Genus)")) +
        geom_bar(stat = "identity", size = 0.1, color = "black") +
        facet_wrap(~Group2, scales = "free_x", nrow = 1) + guides(fill = guide_legend(ncol = 2)) +
        scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
        theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7),
              axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
              axis.text.y = element_text(size = 12)) + xlab("Sample") + ylab("Relative Abundance")

png("images/plot_bar_genus.png", width = 10, height = 8, units = "in", res = 300)
ptaxa6
invisible(dev.off())

#' 

#' 
#' # Alpha diversity
#' 
#' Explore alpha diversity using the `plot_richness` function 
#' 
## ----plot_richness------------------------------------------------------------
# Select alpha-diversity measures
divIdx = c("Observed", "Chao1", "Shannon", "Simpson")

png("images/plot_richness.png", width = 6, height = 5, units = "in", res = 300)
plot_richness(ps1, x = "Group2", measures = divIdx, color = "Group2", nrow = 1) + 
	geom_point(size = 0.8) + theme_bw() + 
	theme(legend.position = "none", 
	      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
	labs(x = "Group", y = "Alpha Diversity Measure")
invisible(dev.off())

#' 

#' 
#' Create a long-format data.frame and plot with `ggplot` function
#' 
## ----plot_richness_boxplot----------------------------------------------------
ad = estimate_richness(ps1, measures = divIdx)
ad = merge(data.frame(sample_data(ps1)), ad, by = "row.names")
ad = ad %>% select(Sample_ID, Group2, all_of(divIdx)) %>% 
	gather(key = "alpha", value = "measure", -c(Sample_ID, Group2))

png("images/plot_richness_boxplot.png", width = 6, height = 5, units = "in", res = 300)
ggplot(ad, aes(Group2, measure, color = Group2)) + 
	geom_boxplot(outlier.shape = NA, size = 0.8, width = 0.8) + 
	geom_quasirandom(size = 0.8, color = "black") + theme_bw() + 
	facet_wrap(~ alpha, scales = "free_y", nrow = 1) + 
	theme(legend.position = "none", 
	      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
	labs(x = "Group", y = "Alpha Diversity Measure")
invisible(dev.off())

#' 

#' 
#' # Beta diversity
#' 
#' Perform various flavors of unifrac measurements using the `GUniFrac` package
#' 
## ----calculate-unifracs, cache = TRUE-----------------------------------------
unifracs = GUniFrac::GUniFrac(otu_table(ps1), phy_tree(ps1), alpha = c(0, 0.5, 1))$unifracs

#' 
#' Perform ordination
#' 
## ----ordinate, cache = TRUE---------------------------------------------------
# Calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm = TRUE){
	exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# Perform variance stabilizing transformation
ps1.vsd = ps1
dds = phyloseq_to_deseq2(ps1, ~ Group2)
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)
dds = estimateDispersions(dds)
abund.vsd = getVarianceStabilizedData(dds)
abund.vsd[abund.vsd < 0] = 0	# set negative values to 0
otu_table(ps1.vsd) = otu_table(abund.vsd, taxa_are_rows = TRUE)

# Create distance objects
dist_un = as.dist(unifracs[, , "d_UW"])       # Unweighted UniFrac
attr(dist_un, "method") = "Unweighted UniFrac"

dist_wu = as.dist(unifracs[, , "d_1"])        # Weighted UniFrac
attr(dist_wu, "method") = "Weighted UniFrac"

dist_vu = as.dist(unifracs[, , "d_VAW"])      # Variance-adjusted-weighted UniFrac
attr(dist_vu, "method") = "Variance-adjusted-weighted UniFrac"

dist_gu = as.dist(unifracs[, , "d_0.5"])      # GUniFrac with alpha 0.5
attr(dist_gu, "method") = "GUniFrac with alpha 0.5"

dist_bc = phyloseq::distance(ps1.vsd, "bray") # Bray-Curtis

# Perform ordination
ord_un = ordinate(ps1, method = "PCoA", distance = dist_un)
ord_wu = ordinate(ps1, method = "PCoA", distance = dist_wu)
ord_vu = ordinate(ps1, method = "PCoA", distance = dist_vu)
ord_gu = ordinate(ps1, method = "PCoA", distance = dist_gu)
set.seed(12345)
ord_bc = ordinate(ps1, method = "NMDS", distance = dist_bc)

#' 
#' Output to PDF
#' 
## ----plot-unifracs, cache = TRUE----------------------------------------------
pdf("images/plot_ordination.pdf", width = 6, height = 5, pointsize = 12)
plot_ordination(ps1, ord_un, type = "samples", color = "Group2", 
		title = "PCoA on unweighted UniFrac") + 
		theme_bw() + coord_fixed(ratio = 1) + 
		stat_ellipse(aes(group = Group2), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_wu, type = "samples", color = "Group2", 
		title = "PCoA on weighted UniFrac") + 
		theme_bw() + coord_fixed(ratio = 1) + 
		stat_ellipse(aes(group = Group2), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_vu, type = "samples", color = "Group2", 
		title = "PCoA on Variance adjusted weighted UniFrac") + 
		theme_bw() + coord_fixed(ratio = 1) + 
		stat_ellipse(aes(group = Group2), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_gu, type = "samples", color = "Group2", 
		title = "PCoA on GUniFrac with alpha 0.5") + 
		theme_bw() + coord_fixed(ratio = 1) + 
		stat_ellipse(aes(group = Group2), type = "t", size = 0.5, linetype = "dashed")

plot_ordination(ps1, ord_bc, type = "samples", color = "Group2", 
		title = "NMDS on Bray-Curtis distance (VST)") + 
		theme_bw() + coord_fixed(ratio = 1) + 
		stat_ellipse(aes(group = Group2), type = "t", size = 0.5, linetype = "dashed")
invisible(dev.off())

#' 
#' Example: PCoA on weighted UniFrac
#' 
## ----plot-weighted-unifrac, echo = FALSE, warning = FALSE, fig.width = 6, fig.height = 5, fig.align = "center", dpi = 100, cache = TRUE----
plot_ordination(ps1, ord_vu, type = "samples", color = "Group2",
                title = "PCoA on weighted UniFrac") + theme_bw() + coord_fixed(ratio = 1) +
                stat_ellipse(aes(group = Group2), type = "t", size = 0.5, linetype = "dashed")

#' 
#' # Multivariate analyses
#' 
#' ## PERMANOVA
#' 
#' Permutational multivariate analysis of variance (PERMANOVA) test for differences between independent groups using the `adonis` function
#' 
## ----adonis-------------------------------------------------------------------
set.seed(12345)
adonis(dist_un ~ Group2, data.frame(sample_data(ps1)), permutations = 1000)

set.seed(12345)
adonis(dist_wu ~ Group2, data.frame(sample_data(ps1)), permutations = 1000)

set.seed(12345)
adonis(dist_vu ~ Group2, data.frame(sample_data(ps1)), permutations = 1000)

set.seed(12345)
adonis(dist_gu ~ Group2, data.frame(sample_data(ps1)), permutations = 1000)

set.seed(12345)
adonis(dist_bc ~ Group2, data.frame(sample_data(ps1)), permutations = 1000)

#' 
#' ## Multivariate dispersion
#' 
#' Test for homogeneity condition among groups using the `betadisper` function
#' 
## ----betadisper---------------------------------------------------------------
disp1 = betadisper(dist_un, group = data.frame(sample_data(ps1))$Group2)
set.seed(12345)
permutest(disp1, permutations = 1000)

disp2 = betadisper(dist_wu, group = data.frame(sample_data(ps1))$Group2) 
set.seed(12345)
permutest(disp2, permutations = 1000)

disp3 = betadisper(dist_vu, group = data.frame(sample_data(ps1))$Group2)
set.seed(12345)
permutest(disp3, permutations = 1000)

disp4 = betadisper(dist_gu, group = data.frame(sample_data(ps1))$Group2)
set.seed(12345)
permutest(disp4, permutations = 1000)

disp5 = betadisper(dist_bc, group = data.frame(sample_data(ps1))$Group2)
set.seed(12345)
permutest(disp5, permutations = 1000)

#' 
#' Example: Dispersion of unweighted UniFrac
#' 
## ----plot-dispersion, fig.width = 6, fig.height = 6, fig.align = "center", dpi = 100----
plot(disp3, hull = FALSE, ellipse = TRUE)

#' 
#' # Save current workspace
#' 
## ----save-image---------------------------------------------------------------
save.image(file = "2_phyloseq_tutorial.RData")

#' 
#' # Session information
#' 
## ----session-info-------------------------------------------------------------
sessionInfo()

