library(GenomicRanges)
library(AnnotationHub)
###'''
#Use the AnnotationHub package to obtain data on "CpG Islands" in the human genome.

# How many islands exists on the autosomes?'''###

ah <- AnnotationHub()

query_results <- query(ah, c("CpG Islands", "Homo sapiens"))

query_results


cpg_islands <- query_results[["AH5086"]]

cpg_islands_data <- cpg_islands

autosomes <- paste0("chr", 1:22)
autosomal_cpg_islands <- cpg_islands_data[seqnames(cpg_islands_data) %in% autosomes]

num_autosomal_cpg_islands <- length(autosomal_cpg_islands)

num_autosomal_cpg_islands


###'''
#How many CpG Islands exists on chromosome 4.
autosome4<- 'chr4'
autosome_4_cpg_islands <-cpg_islands_data[seqnames(cpg_islands_data)==autosome4]
num_autosome_4_cpg_islands<- length(autosome_4_cpg_islands)
num_autosome_4_cpg_islands

###'''
#Obtain the data for the H3K4me3 histone modification for the H1 cell line 
#from Epigenomics Roadmap, using AnnotationHub. Subset these regions to only 
#keep regions mapped to the autosomes (chromosomes 1 to 22).
# How many bases does these regions cover?
hi_query<- query(ah, c("H3K4me3", "EpigenomeRoadMap",'H1 Cells'))
hi_query[2]#choose second subset which has all the keywords mentioned in the question
h3k4me3_data <- hi_query[["AH29884"]]
autosomal_h3k4me3_regions <- h3k4me3_data[seqnames(h3k4me3_data) %in% autosomes]
total_bases_covered <- sum(width(autosomal_h3k4me3_regions))
total_bases_covered



###'''
#Obtain the data for the H3K27me3 histone modification for the H1 cell line 
#from Epigenomics Roadmap, using the AnnotationHub package. Subset these regions 
#to only keep regions mapped to the autosomes. In the return data, each region 
#has an associated "signalValue". 
# What is the mean signalValue across all regions on the standard chromosomes?
Mhi_query<- query(ah, c("H3K4me3"))
h3k4me3_data_sval <- Mhi_query[["AH23273"]]
autosomal_h3k4me3_regions_sval <- h3k4me3_data_sval[seqnames(h3k4me3_data_sval) %in% autosomes]
keeping<- keepStandardChromosomes(autosomal_h3k4me3_regions_sval)
mean(keeping$signalValue)


###'''
#Bivalent regions are bound by both H3K4me3 and H3K27me3.

# Using the regions we have obtained above, how many bases on the standard chromosomes are bivalently marked?
biValent_region = intersect(autosomal_h3k4me3_regions, keeping)
sum(width(biValent_region))

####''''
#We will examine the extent to which bivalent regions overlap CpG Islands.

#how big a fraction (expressed as a number between 0 and 1) of the bivalent
#regions, overlap one or more CpG Islands?
overlaps_bi_cpg <- findOverlaps(biValent_region, cpg_islands)
length(unique(queryHits(overlaps_bi_cpg)))/length(biValent_region)

####''''
#### How big a fraction (expressed as a number between 0 and 1) of the 
#bases which are part of CpG Islands, are also bivalent marked
fraction_cbg = intersect(biValent_region, cpg_islands)
sum(width(fraction_cbg))/sum(width(cpg_islands))

####'''
####' How many bases are bivalently marked within 10kb of CpG Islands?
resized_cpg = resize(cpg_islands, width=20000+width(cpg_islands), fix="center")
intersect_bi_cbg = intersect(biValent_region, resized_cpg)
sum(width(intersect_bi_cbg))

####''''
####' How big a fraction (expressed as a number between 0 and 1) of the human 
####' genome is contained in a CpG Island?
ah_Hgenome = query(ah, c("hg19","Assembly"))
genome = ah_Hgenome[[1]]
genome = dropSeqlevels(genome, "chrX", pruning.mode = "coarse")
genome = dropSeqlevels(genome, "chrY", pruning.mode = "coarse")
genome = keepStandardChromosomes(genome, pruning.mode = "coarse")
genome_size = sum(as.numeric(seqlengths(genome)))
sum(as.numeric(width(cpg_islands)))/genome_size

####'''
####' Compute an odds-ratio for the overlap of bivalent marks with CpG islands.
inOut = matrix(0, ncol = 2, nrow = 2)
colnames(inOut) = c("in", "out")
rownames(inOut) = c("in", "out")
inOut[1,1] = sum(width(intersect(biValent_region, cpg_islands)))
inOut[1,2] = sum(width(setdiff(biValent_region, cpg_islands)))
inOut[2,1] = sum(width(setdiff(cpg_islands, biValent_region)))
inOut[2,2] = genome_size - sum(inOut)

odd_ratio <- inOut[1,1]*inOut[2,2]/(inOut[1,2]*inOut[2,1])
odd_ratio