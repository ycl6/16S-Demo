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
#' **LEfSe:** [GitHub](https://github.com/SegataLab/lefse), [Paper](https://doi.org/10.1186/gb-2011-12-6-r60)
#' 
#' **export2graphlan**: [GitHub](https://github.com/segatalab/export2graphlan), [Documentation](https://github.com/biobakery/graphlan/wiki/export2graphlan---tutorial)
#' 
#' **GraPhlAn**: [GitHub](https://github.com/biobakery/graphlan), [Documentation](https://github.com/biobakery/graphlan/wiki), [Paper](https://doi.org/10.7717/peerj.1029)
#' 
#' **Demo Dataset:** [PRJEB27564](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB27564) from *[Gut Microbiota in Parkinson's Disease: Temporal Stability and Relations to Disease Progression.](https://pubmed.ncbi.nlm.nih.gov/31221587/) EBioMedicine. 2019;44:691-707*
#' 
#' **License:** GPL-3.0
#' 
#' -----
#' 
#' # Workflow (Continues from Part 2)
#' 
#' > **Note:** `export2graphlan` and `GraPhlAn` can only run in Python 2.7.
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
library("grid")
library("tidyverse")
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
#' We are now using new version of LEfSe. Since we are interested in testing the biomarker using the main class labels (control vs. patient), hence the subject line is removed when creating the LEfSe input table.
#' 
## ----prepare-input, cache = TRUE----------------------------------------------
# Add 'subject = TRUE' to include subject in the data.frame
tax1 = lefse_1name_obj(ps1a, sample_data(ps1a)$Group)
lefse1 = lefse_obj(ps1a)
lefse1 = rbind(tax1, lefse1)

tax2 = lefse_1name_obj(ps1b, sample_data(ps1b)$Time)
lefse2 = lefse_obj(ps1b)
lefse2 = rbind(tax2, lefse2)

# Replace unsupported chars with underscore
lefse1$name = gsub(" ","_",lefse1$name)
lefse1$name = gsub("-","_",lefse1$name)
lefse1$name = gsub("/","_",lefse1$name)

lefse2$name = gsub(" ","_",lefse2$name)
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
#' - `expr1`: Patient vs. control at baseline
#' - `expr2`: Patient at baseline vs. patient at followup
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
#' Set the full-paths to `LEfSe`'s python scripts if they are not in `PATH`
#' 
## ----set-lefse-path, eval = FALSE---------------------------------------------
## lefse-format_input = "/path-to-lefse-format_input.py"	# e.g. /home/ngs/lefse/lefse-format_input.py
## run_lefse = "/path-to-run_lefse.py"			# e.g. /home/ngs/lefse/run_lefse.py

#' 

#' 
#' ## Run LEfSe
#' 
#' ### Run LEfSe in terminal/console
#' 
#' ```
#' /path-to-lefse-format_input.py lefse/expr1.lefse_table.txt lefse/expr1.lefse_table.in -c 1 -o 1000000
#' /path-to-lefse-format_input.py lefse/expr2.lefse_table.txt lefse/expr2.lefse_table.in -c 1 -o 1000000
#' 
#' /path-to-run_lefse.py lefse/expr1.lefse_table.in lefse/expr1.lefse_table.res -b 100 -a 1 -l 1
#' /path-to-run_lefse.py lefse/expr2.lefse_table.in lefse/expr2.lefse_table.res -b 100 -a 1 -l 1
#' ```
#' 
#' ### Run LEfSe in R
#' 
## ----run-lefse, cache = TRUE--------------------------------------------------
system2(lefse_format_input, 
	args = c("lefse/expr1.lefse_table.txt", "lefse/expr1.lefse_table.in", 
		 "-c 1", "-o 1000000"))

system2(lefse_format_input, 
	args = c("lefse/expr2.lefse_table.txt", "lefse/expr2.lefse_table.in", 
		 "-c 1", "-o 1000000"))

# set Kruskal-Wallis alpha (-a) to 1 to allow returning of all P-value to perform adjustment later
system2(run_lefse, 
	args = c("lefse/expr1.lefse_table.in", "lefse/expr1.lefse_table.res", 
		 "-b 100", "-a 1", "-l 1"), stdout = TRUE)

system2(run_lefse, 
	args = c("lefse/expr2.lefse_table.in", "lefse/expr2.lefse_table.res", 
		 "-b 100", "-a 1", "-l 1"), stdout = TRUE)

#' 
#' # Multiple testing correction
#' 
#' Perform multiple testing correction on reported P-values
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
#' > **Note:** `export2graphlan` and `GraPhlAn` can only run in Python 2.7. They do not support Python 3.
#' 
#' > **Note:** Update the coloring `lefse/expr1.colors` and `lefse/expr2.colors` if necessary. The colors should be defined in HSV (hue, saturation, value) scale.
#' 
#' Regarding the color file format:
#' 
#' - The spacing between the group name and HSV scales should be **tab**, not space
#' - The spacing between the HSV scales should be **space** not tab
#' - There should be no new line at the end of the file
#' 
#' ```
#' A       150., 100., 53.
#' B       217., 100., 91.
#' C       4., 85., 84.
#' D       39., 100., 100.
#' ```
#' 
#' ## *\*Activate conda environment*
#' 
#' If you have used `conda` to create a Python 2.7-specific environment for `export2graphlan` and `GraPhlAn`, for example an environment called `graphlan`, you can activate it now to run both programs.
#' 
#' > **Note:** Use `conda deactivate` to deactivate from a currently active conda environment if required.
#' 
#' ```
#' conda activate graphlan
#' ```
#' 
#' ## Run `export2graphlan` and `GraPhlAn` in terminal/console
#' 
#' > **Note:** Remember to change the PATHs to each script accordingly.
#' 
#' Patient vs. control at baseline (expr1)
#' 
#' ```
#' # Run export2graphlan
#' /path-to-export2graphlan.py -i lefse/expr1.lefse_table.txt -o lefse/expr1.lefse_table.res.padj \
#' -t lefse/expr1.graphlan_tree.txt -a lefse/expr1.graphlan_annot.txt --external_annotations 2,3,4,5,6 \
#' --fname_row 0 --biomarkers2colors lefse/expr1.colors
#' 
#' # Run graphlan
#' /path-to-graphlan_annotate.py --annot lefse/expr1.graphlan_annot.txt lefse/expr1.graphlan_tree.txt lefse/expr1.graphlan_outtree.txt
#' 
#' # Convert 'lsqb' and 'rsqb' back to square bracket symbols 
#' sed 's/lsqb/[/' lefse/expr1.graphlan_outtree.txt | sed 's/rsqb/]/' > lefse/expr1.graphlan_outtree_sqb.txt
#' sed 's/lsqb/[/' lefse/expr1.lefse_table.res.padj | sed 's/rsqb/]/' > lefse/expr1.lefse_table.res_sqb.padj
#' 
#' # Create cladogram
#' /path-to-graphlan.py --dpi 150 lefse/expr1.graphlan_outtree_sqb.txt lefse/expr1.graphlan.png --external_legends --size 8 --pad 0.2
#' 
#' # Convert to readble output
#' perl /path-to-lefse.pl lefse/expr1.lefse_table.res_sqb.padj lefse/expr1.graphlan_outtree_sqb.txt lefse/expr1.out
#' ```
#' 

#' 
#' Patient at baseline vs. patient at followup (expr2)
#' 
#' ```
#' # Run export2graphlan
#' /path-to-export2graphlan.py -i lefse/expr2.lefse_table.txt -o lefse/expr2.lefse_table.res.padj \
#' -t lefse/expr2.graphlan_tree.txt -a lefse/expr2.graphlan_annot.txt --external_annotations 2,3,4,5,6 \
#' --fname_row 0 --biomarkers2colors lefse/expr2.colors
#' 
#' # Run graphlan
#' /path-to-graphlan_annotate.py --annot lefse/expr2.graphlan_annot.txt lefse/expr2.graphlan_tree.txt lefse/expr2.graphlan_outtree.txt
#' 
#' # Convert 'lsqb' and 'rsqb' back to square bracket symbols
#' sed 's/lsqb/[/' lefse/expr2.graphlan_outtree.txt | sed 's/rsqb/]/' > lefse/expr2.graphlan_outtree_sqb.txt
#' sed 's/lsqb/[/' lefse/expr2.lefse_table.res.padj | sed 's/rsqb/]/' > lefse/expr2.lefse_table.res_sqb.padj
#' 
#' # Create cladogram
#' /path-to-graphlan.py --dpi 150 lefse/expr2.graphlan_outtree_sqb.txt lefse/expr2.graphlan.png --external_legends --size 8 --pad 0.2
#' 
#' # Convert to readble output
#' perl /path-to-lefse.pl lefse/expr2.lefse_table.res_sqb.padj lefse/expr2.graphlan_outtree_sqb.txt lefse/expr2.out
#' ```
#' 

#' 

#' 
#' > **Note:** Use `conda deactivate` to deactivate the `graphlan` environment if required.
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
png("lefse/expr1.lefse_table.png", width = 6, height = 5, units = "in", res = 300)
ggplot(res1, aes(taxon, LDA, fill = Group)) + geom_bar(stat = "identity", width = 0.7, size = 0.5) + 
	coord_flip() + theme_bw() + facet_grid(rows = vars(Group), scales = "free_y", space = "free_y") +
	scale_fill_manual(values = c("control" = "blue", "patient" = "red")) +
	theme(legend.position = "none", 
	      axis.text.y = element_text(face = "bold", size = 8), 
	      strip.placement = "outside",
	      strip.text.y = element_text(face = "bold", angle = 0)) +
	labs(title = "LEfSe Analysis (Baseline)", x = "Taxon", y = "LDA")
invisible(dev.off())

#' 

#' 
#' Patient at baseline vs. patient at followup (expr2)
#' 
## ----plot-bar-expr2, cache = TRUE---------------------------------------------
png("lefse/expr2.lefse_table.png", width = 6, height = 3, units = "in", res = 300)
ggplot(res2, aes(taxon, LDA, fill = Group)) + geom_bar(stat = "identity", width = 0.7, size = 0.5) + 
	coord_flip() + theme_bw() + facet_grid(rows = vars(Group), scales = "free_y", space = "free_y") +
        scale_fill_manual(values = c("baseline" = "blue", "followup" = "red")) +
        theme(legend.position = "none", 
              axis.text.y = element_text(face = "bold", size = 8), 
	      strip.placement = "outside",
              strip.text.y = element_text(face = "bold", angle = 0)) +
        labs(title = "LEfSe Analysis (followup)", x = "Taxon", y = "LDA")
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

