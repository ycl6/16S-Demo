#' ---
#' title: "16S rDNA amplicon sequencing analysis using R (Part 3: LEfSe & GraPhlAn)"
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
#' **LEfSe:** [Bitbucket](https://bitbucket.org/nsegata/lefse/downloads/), [Paper](https://doi.org/10.1186/gb-2011-12-6-r60)
#' 
#' **export2graphlan**: [GitHub](https://github.com/segatalab/export2graphlan), [export2graphlan](https://bitbucket.org/nsegata/graphlan/wiki/export2graphlan%20-%20tutorial)
#' 
#' **GraPhlAn**: [GitHub](https://github.com/SegataLab/graphlan), [Documentation](https://bitbucket.org/nsegata/graphlan/wiki/Home), [Paper](https://doi.org/10.7717/peerj.1029)
#' 
#' **Demo Dataset:** [PRJEB27564](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB27564) from *[Gut Microbiota in Parkinson's Disease: Temporal Stability and Relations to Disease Progression.](https://pubmed.ncbi.nlm.nih.gov/31221587/) EBioMedicine. 2019;44:691-707*
#' 
#' **License:** GPL-3.0
#' 
#' -----
#' 
#' # Workflow (Continues from Part 2)
#' 
#' > **Note:** Requires Python 2
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
library("ggplot2")
library("dplyr")
library("grid")
Sys.setlocale("LC_COLLATE", "C")

#' 
#' Set the full-path to additional scripts
#' 
## ----set-lefse-Rscript-path, eval = FALSE-------------------------------------
## source("/path-to-script/lefse.R", local = TRUE)

#' 

#' 
#' # Load saved workspace from Part 2
#' 
## ----load-image---------------------------------------------------------------
load("2_phyloseq_tutorial.RData")

#' 
#' Create a folder `lefse`
#' 
## ----dir-create-lefse---------------------------------------------------------
if(!dir.exists("lefse")) { dir.create("lefse") }

#' 
#' # Prepare input data
#' 
#' ## Subset samples
#' 
## ----subset-samples, cache = TRUE---------------------------------------------
# Patient vs. control at baseline
ps1a = prune_samples(sample_names(ps1)[grep("B$", sample_names(ps1))], ps1)
ps1a = prune_taxa(taxa_sums(ps1a) > 0, ps1a)

# Patient followup vs. patient baseline
ps1b = prune_samples(sample_names(ps1)[grep("^P", sample_names(ps1))], ps1)
ps1b = prune_taxa(taxa_sums(ps1b) > 0, ps1b)

#' 
#' ## Prepare LEfSe input
#' 
## ----prepare-input, cache = TRUE----------------------------------------------

# Patient vs. control at baseline
ps1a = prune_samples(sample_names(ps1)[grep("B$", sample_names(ps1))], ps1)
ps1a = prune_taxa(taxa_sums(ps1a) > 0, ps1a)

# Patient followup vs. patient baseline
ps1b = prune_samples(sample_names(ps1)[grep("^P", sample_names(ps1))], ps1)
ps1b = prune_taxa(taxa_sums(ps1b) > 0, ps1b)

tax1 = lefse_1name_obj(ps1a, sample_data(ps1a)$Group)
lefse1 = lefse_obj(ps1a)
lefse1 = rbind(tax1, lefse1)

tax2 = lefse_1name_obj(ps1b, sample_data(ps1b)$Time)
lefse2 = lefse_obj(ps1b)
lefse2 = rbind(tax2, lefse2)

# Replace unsupported chars
lefse1$name = gsub("-","_",lefse1$name)
lefse1$name = gsub("/","_",lefse1$name)

lefse2$name = gsub("-","_",lefse2$name)
lefse2$name = gsub("/","_",lefse2$name) 

#' 
#' > For **Silva v132** users, apply below to fix identifical phylum & class names (Actinobacteria and Deferribacteres)
#' 
## ----fix-taxa-silva231, eval = FALSE------------------------------------------
## lefse1 = fix_taxa_silva132(lefse1)
## lefse2 = fix_taxa_silva132(lefse2)

#' 
#' ## Output to file
#' 
## ----output-lefse-table-------------------------------------------------------
write.table(lefse1, file = "lefse/expr1.lefse_table.txt", sep = "\t", quote = F, 
	    row.names = F, col.names = F)
write.table(lefse2, file = "lefse/expr2.lefse_table.txt", sep = "\t", quote = F, 
	    row.names = F, col.names = F)

#' 
#' # Perform LEfSe analysis
#' 
#' We use linear discriminant analysis (LDA) effect size (LEfSe) to determines the features (in this case the clades at each taxonomic rank) most likely to explain the differences between Parkinson's patients and control subjects
#' 
#' ## Set the full-path to LEfSe folder
#' 
## ----set-nsegata-lefse-path, eval = FALSE-------------------------------------
## nsegata_lefse = "/path-to-script/nsegata-lefse-9adc3a62460e/"

#' 

#' 
#' ## Run LEfSe
#' 
## ----run-lefse, cache = TRUE--------------------------------------------------
system2(paste0(nsegata_lefse, "format_input.py"), 
	args = c("lefse/expr1.lefse_table.txt", "lefse/expr1.lefse_table.in", 
		 "-c 1", "-u 2", "-o 1000000"))

system2(paste0(nsegata_lefse, "format_input.py"), 
	args = c("lefse/expr2.lefse_table.txt", "lefse/expr2.lefse_table.in", 
		 "-c 1", "-u 2", "-o 1000000"))

# set KW alpha to 1 to allow returning of all P-value to perform adjustment later
system2(paste0(nsegata_lefse, "run_lefse.py"), 
	args = c("lefse/expr1.lefse_table.in", "lefse/expr1.lefse_table.res", 
		 "-b 100", "-a 1", "-l 2"), stdout = TRUE)

system2(paste0(nsegata_lefse, "run_lefse.py"), 
	args = c("lefse/expr2.lefse_table.in", "lefse/expr2.lefse_table.res", 
		 "-b 100", "-a 1", "-l 2"), stdout = TRUE)

#' 
#' # Multiple testing correction
#' 
#' Perform multiple testing correction on KW P-values
#' 
## ----p-correction, warning = FALSE--------------------------------------------
# set fdr threshold at 0.1
q = 0.1

expr1 = data.frame(fread("lefse/expr1.lefse_table.res"))
expr1 = pcorrection(expr1, q)

expr2 = data.frame(fread("lefse/expr2.lefse_table.res"))
expr2 = pcorrection(expr2, q)

# You can use table() to check number of taxa pass the threshold
table(expr1$V3 != "")

table(expr2$V3 != "")

# Write new result file with the P-value column replaced by the FDR
write.table(expr1[,c(1:4,6)], file = "lefse/expr1.lefse_table.res.padj", sep = "\t", quote = F, 
	    row.names = F, col.names = F)
write.table(expr2[,c(1:4,6)], file = "lefse/expr2.lefse_table.res.padj", sep = "\t", quote = F, 
	    row.names = F, col.names = F)

#' 
#' # Plot cladogram
#' 
#' ## Set paths
#' 
#' Set the full-paths to `export2graphlan`, `GraPhlAn` and `lefse.pl` if they are not in PATH
#' 
## ----set-graphlan-path, eval = FALSE------------------------------------------
## export2graphlan = "/path-to-script/export2graphlan.py"
## graphlan_annotat = "/path-to-script/graphlan_annotate.py"
## graphlan = "/path-to-script/graphlan.py"
## graphlan_parser = "/path-to-script/lefse.pl"

#' 

#' 
#' ## Run `export2graphlan` and `GraPhlAn`
#' 
#' > **Note:** Update the coloring `lefse/expr1.colors` and `lefse/expr2.colors` if necessary. The colors should be defined in HSV (hue, saturation, value) scale
#' 
#' Patient vs. control at baseline (expr1)
#' 
## ----plot-cladogram-1, cache = TRUE-------------------------------------------
system2(export2graphlan, 
	args = c("-i lefse/expr1.lefse_table.txt", "-o lefse/expr1.lefse_table.res.padj", 
		 "-t lefse/expr1.graphlan_tree.txt", "-a lefse/expr1.graphlan_annot.txt", 
		 "--external_annotations 2,3,4,5,6", "--fname_row 0", "--skip_rows 1", 
		 "--biomarkers2colors lefse/expr1.colors"))

system2(graphlan_annotat, 
	args = c("--annot lefse/expr1.graphlan_annot.txt", 
		 "lefse/expr1.graphlan_tree.txt", 
		 "lefse/expr1.graphlan_outtree.txt"))

system2(graphlan, 
	args = c("--dpi 150", "lefse/expr1.graphlan_outtree.txt", "lefse/expr1.graphlan.png", 
		 "--external_legends", "--size 8", "--pad 0.2"))

#' 

#' 
#' Patient at baseline vs. patient at followup (expr2)
#' 
## ----plot-cladogram-2, cache = TRUE-------------------------------------------
system2(export2graphlan, 
        args = c("-i lefse/expr2.lefse_table.txt", "-o lefse/expr2.lefse_table.res.padj", 
                 "-t lefse/expr2.graphlan_tree.txt", "-a lefse/expr2.graphlan_annot.txt", 
                 "--external_annotations 2,3,4,5,6", "--fname_row 0", "--skip_rows 1", 
                 "--biomarkers2colors lefse/expr2.colors"))

system2(graphlan_annotat, 
        args = c("--annot lefse/expr2.graphlan_annot.txt", 
                 "lefse/expr2.graphlan_tree.txt",
                 "lefse/expr2.graphlan_outtree.txt"))

system2(graphlan, 
        args = c("--dpi 150", "lefse/expr2.graphlan_outtree.txt", "lefse/expr2.graphlan.png", 
                 "--external_legends", "--size 8", "--pad 0.2"))

#' 

#' 
#' ## Convert to readble output
#' 
## ----parse-graphlan, cache = TRUE---------------------------------------------
system2("perl", args = c(graphlan_parser, "lefse/expr1.lefse_table.res.padj", 
			 "lefse/expr1.graphlan_outtree.txt", "lefse/expr1.out"))

system2("perl", args = c(graphlan_parser, "lefse/expr2.lefse_table.res.padj", 
			 "lefse/expr2.graphlan_outtree.txt", "lefse/expr2.out"))

#' 
#' # Plot results
#' 
#' ## Load parsed results
#' 
#' Patient vs. control at baseline (expr1)
#' 
## ----load-res-expr1, cache = TRUE---------------------------------------------
res1 = data.frame(fread("lefse/expr1.out"))
res1 = res1[order(res1$order),]
names(res1)[c(5:7)] = c("Group","LDA","FDR")
res1$taxon = paste0(res1$taxon, "(",res1$rank, ";", " ",res1$order, ")")
res1$taxon = as.factor(res1$taxon)
res1$taxon = factor(res1$taxon, levels = unique(res1$taxon[order(res1$order, decreasing = T)]))
res1$Group = as.factor(res1$Group)
levels(res1$taxon) = gsub("Escherichia_Shigella","Escherichia/Shigella",levels(res1$taxon))
levels(res1$taxon) = gsub("_UCG_","_UCG-",levels(res1$taxon))

#' 

#' 
#' Patient at baseline vs. patient at followup (expr2)
#' 
## ----load-res-expr2, cache = TRUE---------------------------------------------
res2 = data.frame(fread("lefse/expr2.out"))
res2 = res2[order(res2$order),]
names(res2)[c(5:7)] = c("Group","LDA","FDR")
res2$taxon = paste0(res2$taxon, "(",res2$rank, ";", " ",res2$order, ")")
res2$taxon = as.factor(res2$taxon)
res2$taxon = factor(res2$taxon, levels = unique(res2$taxon[order(res2$order, decreasing = T)]))
res2$Group = as.factor(res2$Group)
levels(res2$taxon) = gsub("Escherichia_Shigella","Escherichia/Shigella",levels(res2$taxon))
levels(res2$taxon) = gsub("_UCG_","_UCG-",levels(res2$taxon))

#' 

#' 
#' ## Plot bar
#' 
#' Patient vs. control at baseline (expr1)
#' 
## ----plot-bar-expr1, cache = TRUE---------------------------------------------
plot_lefse = ggplot(res1, aes(taxon, LDA, fill = Group)) + 
	geom_bar(stat = "identity", width = 0.7, size = 0.5) + coord_flip() + 
	theme_bw() + facet_wrap(~ Group, ncol = 1, scales = "free_y") + 
	scale_fill_manual(values = c("control" = "blue", "patient" = "red")) +
	theme(legend.position = "none", 
	      axis.text.x = element_text(size = 12), 
	      axis.text.y = element_text(face = "bold", size = 8), 
	      strip.text.x = element_text(face = "bold", size = 12)) + 
	labs(title = "LEfSe Analysis (Baseline)", x = "Taxon", y = "LDA")

groupN = res1 %>% group_by(Group) %>% summarise(count = length(unique(taxon)))
gt = ggplotGrob(plot_lefse)
panelI.1 = gt$layout$t[grepl("panel", gt$layout$name)]
gt$heights[panelI.1] = unit(groupN$count, "null")
invisible(dev.off())

png("lefse/expr1.lefse_table.png", width = 6, height = 6, units = "in", res = 300)
grid.draw(gt)
invisible(dev.off())

#' 

#' 
#' Patient at baseline vs. patient at followup (expr2)
#' 
## ----plot-bar-expr2, cache = TRUE---------------------------------------------
plot_lefse = ggplot(res2, aes(taxon, LDA, fill = Group)) + 
        geom_bar(stat = "identity", width = 0.7, size = 0.5) + coord_flip() + 
        theme_bw() + facet_wrap(~ Group, ncol = 1, scales = "free_y") +  
        scale_fill_manual(values = c("baseline" = "blue", "followup" = "red")) +
        theme(legend.position = "none", 
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(face = "bold", size = 8), 
              strip.text.x = element_text(face = "bold", size = 12)) + 
        labs(title = "LEfSe Analysis (followup)", x = "Taxon", y = "LDA")

groupN = res2 %>% group_by(Group) %>% summarise(count = length(unique(taxon)))
gt = ggplotGrob(plot_lefse)
panelI.1 = gt$layout$t[grepl("panel", gt$layout$name)]
gt$heights[panelI.1] = unit(groupN$count, "null")
invisible(dev.off())

png("lefse/expr2.lefse_table.png", width = 6, height = 3, units = "in", res = 300)
grid.draw(gt)
invisible(dev.off())

#' 

#' 
#' # Save current workspace
#' 
## ----save-image---------------------------------------------------------------
save.image(file = "3_lefse_tutorial.RData")

#' 
#' # Session information
#' 
## ----session-info-------------------------------------------------------------
sessionInfo()

