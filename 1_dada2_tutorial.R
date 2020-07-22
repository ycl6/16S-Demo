#' ---
#' title: "16S rDNA amplicon sequencing analysis using R (Part 1: cutadapt & DADA2)"
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
#' **cutadapt:** [GitHub](https://github.com/marcelm/cutadapt), [Documentation](https://cutadapt.readthedocs.io/en/stable/guide.html), [Paper](https://doi.org/10.14806/ej.17.1.200)
#' 
#' **DADA2:** [GitHub](https://github.com/benjjneb/dada2), [Documentation](https://benjjneb.github.io/dada2/index.html), [Paper](https://doi.org/10.1038/nmeth.3869)
#' 
#' **Demo Dataset:** [PRJEB27564](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB27564) from *[Gut Microbiota in Parkinson's Disease: Temporal Stability and Relations to Disease Progression.](https://pubmed.ncbi.nlm.nih.gov/31221587/) EBioMedicine. 2019;44:691-707*
#' 
#' **License:** GPL-3.0
#' 
#' -----
#' 
#' # Primer trimming using `cutadapt`
#' 
#' You can use the included Perl script `run_trimming.pl` to perform primer trimming, where the raw fastq files are placed in the `raw` folder
#' 
#' ```
#' cd /ngs/16S-Demo
#' 
#' # Usage: perl run_trimming.pl project_folder fastq_folder forward_primer_seq reverse_primer_seq
#' perl run_trimming.pl PRJEB27564 raw CCTACGGGNGGCWGCAG GACTACHVGGGTATCTAATCC
#' ```
#' 
#' Or follow the instruction below to run `cutadapt` within R
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
library("dada2")
library("data.table")
library("phyloseq")
library("ggplot2")

#' 
#' Change the `fastq` path to correspond to the location of the raw fastq files
#' 
## ----set-path-----------------------------------------------------------------
fastq = "raw"           # raw fastq files
trimmed = "trimmed"     # cutadapt trimmed fastq files
filt = "filt"           # dada2 trimmed fastq files
outfiles = "outfiles"   # output files
images = "images"       # output images

if(!dir.exists(filt)) dir.create(filt)
if(!dir.exists(outfiles)) dir.create(outfiles)
if(!dir.exists(images)) dir.create(images)

#' 
## -----------------------------------------------------------------------------
head(list.files(fastq),15)

#' 
#' # *\*Primer trimming in R*
#' 
#' The trimming program `cutadapt` is called using `system2` function to perform trimming.
#' 
#' > **Note:** Skip this step if trimmed has been performed
#' 
## ----trimming, eval = FALSE---------------------------------------------------
## fns = sort(list.files(fastq, full.names = TRUE))
## fnFs = fns[grep("1.fastq.gz", fns)]
## fnRs = fns[grep("2.fastq.gz", fns)]
## 
## if(!dir.exists(trimmed)) dir.create(trimmed)
## 
## fnFs.cut = file.path(trimmed, basename(fnFs))
## fnRs.cut = file.path(trimmed, basename(fnRs))
## log.cut = gsub(".1.fastq.gz", ".log", fnFs.cut)
## sample.names = gsub(".1.fastq.gz", "", basename(fnFs.cut))
## 
## # Define the primer set used to perform PCR
## FWD = "CCTACGGGNGGCWGCAG"       # 341F
## REV = "GACTACHVGGGTATCTAATCC"   # 785R
## 
## # Get reverse complement DNA sequences
## FWD.RC = dada2::rc(FWD)
## REV.RC = dada2::rc(REV)
## 
## # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
## R1.flags = paste("-g", FWD, "-a", REV.RC)
## # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
## R2.flags = paste("-G", REV, "-A", FWD.RC)
## 
## # Run cutadapt to remove primers
## # Note to change the PATH to cutadapt accordingly
## cutadapt = "/path-to-bin/cutadapt"
## 
## for(i in seq_along(fnFs)) {
## 	print(paste0("[", i ,"/", length(sample.names), "] ", sample.names[i]))
## 	
## 	system2(cutadapt,
## 		stdout = log.cut[i], stderr = log.cut[i], # log file
## 		args = c(R1.flags, R2.flags,
## 		"-n 2",                   # -n 2 required to remove FWD and REV from reads
## 		"--match-read-wildcards", # enable IUPAC nucleotide codes (wildcard characters)
## 		"--length 300",           # Truncate reads to 300 bp
## 		"-m 150",                 # discard reads shorter than LEN (avoid length zero sequences)
## 		"--overlap 10",           # min overlap between read and adapter for an adapter to be found
## 		"-j 0",                   # auto-detection of CPU cores, only available on Python 3
## 		"-o", fnFs.cut[i], "-p", fnRs.cut[i], # trimmed files
## 		fnFs[i], fnRs[i])         # input files
## 	)
## }

#' 
## -----------------------------------------------------------------------------
head(list.files(trimmed),15)

#' 
#' # Build trimmed file lists
#' 
## ----build-trimmed-file-lists-------------------------------------------------
fns = sort(list.files(trimmed, full.names = TRUE))
fnFs = fns[grep("1.fastq.gz", fns)]     # Update the grep pattern when necessary
fnRs = fns[grep("2.fastq.gz", fns)]     # Update the grep pattern when necessary
sample.names = gsub(".1.fastq.gz", "", basename(fnFs))  # Update the gsub pattern when necessary

# Check objects
fnFs
fnRs
sample.names

#' 
#' # Inspect read quality profiles
#' 
#' We use the `plotQualityProfile` function to plot the quality profiles of the trimmed fastq files
#' 
## ----plotQualityProfile, message = FALSE, cache = TRUE------------------------
# Plot quality profile of fastq files
ii = 1:length(sample.names)
pdf(paste0(images, "/plotQualityProfile.pdf"), width = 8, height = 8, pointsize = 12)
for(i in ii) {
	message(paste0("[", i ,"/", length(sample.names), "] ", sample.names[i]))
	print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd"))
	print(plotQualityProfile(fnRs[i]) + ggtitle("Rev"))
}
invisible(dev.off())

#' 
#' We review the quality profiles, the quality distribution of forward reads appears to drop after position 260 and position 200 for reverse reads. We will use these positions to truncate forward and reverse reads respectively in the next step.
#' 

#' 
#' # Filter and trim
#' 
#' Remember to review **plotQualityProfile.pdf** to select the best paramters for `truncLen` argument
#' 
## ----filterAndTrim, message = FALSE, cache = TRUE-----------------------------
# Set paths to the dada2-filterd files
filtFs = file.path(filt, basename(fnFs))
filtRs = file.path(filt, basename(fnRs))

# Perform filtering and trimming
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs,
        # Need to keep paramters consistent between runs of the same study
        truncLen = c(260,200), minLen = 200, maxN = 0, truncQ = 2, maxEE = c(2,5),
        rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE)

out = as.data.frame(out)
rownames(out) = sample.names

#' 
## -----------------------------------------------------------------------------
head(out, 10)

#' 
#' # Learn the Error Rates
#' 
#' Use the `learnErrors` function to perform dereplication and learn the error rates. The `derepFastq` function used in past workflow has been intergrated into `learnErrors` function
#' 
## ----learnErrors, cache = TRUE------------------------------------------------
errF = learnErrors(filtFs, multithread = TRUE)

errR = learnErrors(filtRs, multithread = TRUE)

# Visualize the estimated error rates
pdf(paste0(images, "/plotErrors.pdf"), width = 10, height = 10, pointsize = 12)
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
invisible(dev.off())

#' 
#' Forward reads
#' 

#' 
#' Reverse reads
#' 

#' 
#' # Sample Inference
#' 
#' We now apply the core sample inference algorithm to the filtered and trimmed sequence data
#' 
#' > **Note:** By default, the `dada` function processes each sample independently (i.e. `pool = FALSE`), if your samples are from an extremely diverse community (e.g. soil), pooling (i.e. `pool = TRUE`) or pseudo-pooling (recommended; i.e. `pool = pseudo`) might help in identifying the rare ASVs in each sample
#' 
## ----dada-inference, results = "hide", cache = TRUE---------------------------
dadaFs = dada(filtFs, err = errF, pool = FALSE, multithread = TRUE)
dadaRs = dada(filtRs, err = errR, pool = FALSE, multithread = TRUE)

#' 
#' # Merge paired reads
#' 
#' We now merge the forward and reverse reads together to obtain the full denoised sequences
#' 
## ----mergePairs, message = FALSE, cache = TRUE--------------------------------
mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

#' 
#' # Construct sequence table
#' 
#' We now construct an amplicon sequence variant (ASV) table
#' 
## ----makeSequenceTable--------------------------------------------------------
seqtab = makeSequenceTable(mergers)

#' 
#' View the length frequency distribution with `table`
#' 
## -----------------------------------------------------------------------------
table(nchar(getSequences(seqtab)))

#' 
#' Save sequence table
#' 
## ----save-rds-----------------------------------------------------------------
saveRDS(seqtab, "seqtab.rds")

# Or as an example, save only the first 5 samples
# saveRDS(seqtab[c(1:5),], "seqtab.rds")

#' 
#' # *\*Merge multiple runs*
#' 
#' This section elaborates a little about merging results from multiple sequencing run. You can merge sequence tables from multiple runs belonging to the same experiment or project before moving on to chimera removal and taxonomy assignment. In the example snippet below, we use `mergeSequenceTables` to combine 2 sequence tables.
#' 
#' > **Note:** Skip this step if you do not require merging multiple sequence tables
#' 
## ----mergeSequenceTables, eval = FALSE----------------------------------------
## # Load sequence tables from multiple runs
## seqtab1 = readRDS("/some-path/run1/seqtab.rds")
## seqtab2 = readRDS("/some-path/run2/seqtab.rds")
## 
## # Check the dimensions and sample names
## dim(seqtab1)
## rownames(seqtab1)
## dim(seqtab2)
## rownames(seqtab2)
## 
## # Merge sequence table and remove chimeras
## seqtab = mergeSequenceTables(seqtab1, seqtab2)
## 
## # To sum values of same sample from multiple sequence table (i.e. when a sample was re-sequenced due to low depth)
## # seqtab = mergeSequenceTables(seqtab1, seqtab2, repeats = "sum")
## 
## # Check the dimension of the merged sequence table
## dim(seqtab)

#' 
#' After obtaining the merged sequence table, you can continue with the below data processing steps
#' 
#' # Remove chimeras
#' 
## ----removeBimeraDenovo, cache = TRUE-----------------------------------------
seqtab.nochim = removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, 
				   verbose = TRUE)

# Update the same names to exclude "fastq.gz" in name
rownames(seqtab.nochim) = sample.names

#' 
#' View the dimensions of `seqtab` and `seqtab.nochim`
#' 
## -----------------------------------------------------------------------------
dim(seqtab)
dim(seqtab.nochim)

#' 
#' Print the proportion of non-chimeras in merged sequence reads
#' 
## -----------------------------------------------------------------------------
sum(seqtab.nochim)/sum(seqtab)

#' 
#' # Track reads through the pipeline
#' 
## ----create-track, cache = TRUE-----------------------------------------------
getN <- function(x) sum(getUniques(x))

track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
	      sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) = c("Trimmed", "Filtered", "denoisedF", "denoisedR", "merged", "nonchim")
track = cbind(data.frame(SampleID = sample.names), track)

write.table(track, file = paste0(outfiles, "/track.txt"), sep = "\t", quote = F, 
	    row.names = F, col.names = T)

#' 

#' 
#' # Assign taxonomy
#' 
#' Set the full-path to the Silva and NCBI database
#' 
## ----set-db-path, eval = FALSE------------------------------------------------
## dbpath = "/path-to-db/"
## ref1 = paste0(dbpath, "silva_nr_v138_train_set.fa.gz")
## ref2 = paste0(dbpath, "silva_species_assignment_v138.fa.gz")
## ref3 = paste0(dbpath, "16SMicrobial.fa.gz")
## 

#' 

#' 
#' Use the `assignTaxonomy` function to classifies sequences against **SILVA** reference training dataset `ref1`, and use the `assignSpecies` function to perform taxonomic assignment to the species level by exact matching against **SILVA** `ref2` and **NCBI** reference datasets `ref3`
#' 
## ----assignTaxonomy, cache = TRUE---------------------------------------------
taxtab = assignTaxonomy(seqtab.nochim, refFasta = ref1, minBoot = 80, 
			tryRC = TRUE, outputBootstraps = TRUE, verbose = TRUE, multithread = TRUE)

spec_silva = assignSpecies(getSequences(seqtab.nochim), ref2, allowMultiple = FALSE, 
			   tryRC = TRUE, verbose = TRUE)

spec_ncbi = assignSpecies(getSequences(seqtab.nochim), ref3, allowMultiple = FALSE, 
			  tryRC = TRUE, verbose = TRUE)

#' 
#' Combine species-level taxonomic assignment from 2 reference sources
#' 
## ----combine-assignment-------------------------------------------------------
SVformat = paste("%0",nchar(as.character(ncol(seqtab.nochim))),"d", sep = "")
svid = paste0("SV", sprintf(SVformat, seq(ncol(seqtab.nochim))))

s_silva = as.data.frame(spec_silva, stringsAsFactors = FALSE)
rownames(s_silva) = svid

s_ncbi = as.data.frame(spec_ncbi, stringsAsFactors = FALSE)
rownames(s_ncbi) = svid
s_ncbi$Genus = gsub("\\[|\\]", "", s_ncbi$Genus)

s_merged = cbind(s_ncbi, s_silva)
colnames(s_merged) = c("nGenus","nSpecies","sGenus","sSpecies")
s_merged1 = s_merged[!is.na(s_merged$nSpecies),]
colnames(s_merged1)[1:2] = c("Genus","Species")
s_merged2 = s_merged[is.na(s_merged$nSpecies) & !is.na(s_merged$sSpecies),]
colnames(s_merged2)[3:4] = c("Genus","Species")
s_merged3 = s_merged[is.na(s_merged$nSpecies) & is.na(s_merged$sSpecies),]
colnames(s_merged3)[3:4] = c("Genus","Species")

s_final = rbind(s_merged1[,c("Genus","Species")], s_merged2[,c("Genus","Species")], 
		s_merged3[,c("Genus","Species")])
s_final = s_final[order(row.names(s_final)),]

s_final = as.matrix(s_final)

if("Genus" %in% colnames(taxtab$tax)) {
	gcol = which(colnames(taxtab$tax) == "Genus")
} else {
	gcol = ncol(taxtab$tax)
}

matchGenera <- function(gen.tax, gen.binom, split.glyph = "/") {
	if(is.na(gen.tax) || is.na(gen.binom)) { return(FALSE) }
	if((gen.tax == gen.binom) || 
	   grepl(paste0("^", gen.binom, "[ _", split.glyph, "]"), gen.tax) || 
	   grepl(paste0(split.glyph, gen.binom, "$"), gen.tax)) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}

gen.match = mapply(matchGenera, taxtab$tax[,gcol], s_final[,1])
taxtab$tax = cbind(taxtab$tax, s_final[,2])
colnames(taxtab$tax)[ncol(taxtab$tax)] = "Species"

print(paste(sum(!is.na(s_final[,2])), "out of", 
	    nrow(s_final), "were assigned to the species level."))

taxtab$tax[!gen.match,"Species"] = NA
print(paste("Of which", sum(!is.na(taxtab$tax[,"Species"])), 
	    "had genera consistent with the input table."))

#' 
#' # Multiple sequence alignment
#' 
#' Prepare a `data.frame` `df` from `seqtab.nochim` and `taxtab`
#' 
## ----prepare-df---------------------------------------------------------------
df = data.frame(sequence = colnames(seqtab.nochim), abundance = colSums(seqtab.nochim), 
		stringsAsFactors = FALSE)
df$id = svid
df = merge(df, as.data.frame(taxtab), by = "row.names")
rownames(df) = df$id
df = df[order(df$id),2:ncol(df)]

#' 
#' Performs alignment of multiple unaligned sequences
#' 
## ----alignseqs, cache = TRUE--------------------------------------------------
alignment = DECIPHER::AlignSeqs(Biostrings::DNAStringSet(setNames(df$sequence, df$id)), anchor = NA)

#' 
#' Export alignment
#' 
## ----export-alignment, cache = TRUE-------------------------------------------
phang.align = phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
phangorn::write.phyDat(phang.align, file = "alignment.fasta", format = "fasta")
phangorn::write.phyDat(phang.align, file = "alignment.aln", format = "phylip")

#' 
#' # Construct phylogenetic tree
#' 
#' Set the full-path to the RAxML and RAxML-NG
#' 
## ----set-raxml-path, eval = FALSE---------------------------------------------
## raxml = "/path-to-raxml"
## raxmlng = "/path-to-raxml-ng"

#' 

#' 
#' Run RAxML and RAxML-NG
#' 
## ----run-raxml, results = "hide", cache = TRUE--------------------------------
system2(raxml, args = c("-T 2", "-f E", "-p 1234", "-x 5678", "-m GTRCAT", "-N 1", 
			"-s alignment.aln", "-n raxml_tree_GTRCAT"))

system2(raxmlng, args = c("--evaluate", "--force", "--seed 1234", "--log progress", "--threads 2", 
			  "--msa alignment.fasta",  "--model GTR+G", "--brlen scaled", 
			  "--tree RAxML_fastTree.raxml_tree_GTRCAT", "--prefix GTRCAT"))

#' 
#' Import tree using the `read_tree` function from `phyloseq`
#' 
## ----import-tree--------------------------------------------------------------
raxml_tree = read_tree("GTRCAT.raxml.bestTree")

#' 
#' # Load sample information
#' 
## ----load-sample-data---------------------------------------------------------
samdf = data.frame(fread("sample.meta", colClasses = "character"))

rownames(samdf) = samdf$Sample_ID
samdf$Sample_ID = as.factor(samdf$Sample_ID)
samdf$Sample_ID = factor(samdf$Sample_ID, levels = c(sort(levels(samdf$Sample_ID), decreasing = F)))
samdf$Subject = as.factor(samdf$Subject)
samdf$Group2 = paste(samdf$Group, samdf$Time, sep = "_")
samdf$Group = as.factor(samdf$Group)
samdf$Group2 = as.factor(samdf$Group2)
samdf$Time = as.factor(samdf$Time)

head(samdf)

#' 
#' # Handoff to `phyloseq`
#' 
#' Prepare `new_seqtab` and `tax` data.frames containing abundance and taxonomy information respectively
#' 
## ----create-tax---------------------------------------------------------------
new_seqtab = seqtab.nochim
colnames(new_seqtab) = df[match(colnames(new_seqtab), df$sequence),]$id

# Update rownames in new_seqtab from Run_ID to Sample_ID
rownames(new_seqtab) = as.character(samdf[match(rownames(new_seqtab), samdf$Run_ID),]$Sample_ID)

new_taxtab = taxtab
rownames(new_taxtab$tax) = df[match(rownames(new_taxtab$tax), df$sequence),]$id

tax = as.data.frame(new_taxtab$tax)
tax$Family = as.character(tax$Family)
tax$Genus = as.character(tax$Genus)

#' 
#' Construct a `phyloseq` object
#' 
## ----create-ps-obj------------------------------------------------------------
ps = phyloseq(tax_table(as.matrix(tax)), 
	      sample_data(samdf), 
	      otu_table(new_seqtab, taxa_are_rows = FALSE), 
	      phy_tree(raxml_tree))

ps

#' 
#' # Save current workspace
#' 
## ----save-image---------------------------------------------------------------
save.image(file = "1_dada2_tutorial.RData")

#' 
#' # Session information
#' 
## ----session-info-------------------------------------------------------------
sessionInfo()

