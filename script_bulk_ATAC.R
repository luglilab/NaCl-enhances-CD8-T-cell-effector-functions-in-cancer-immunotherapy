# Code Figure: 5J, 5K, 5I

###################
###################
# Bulk NaCl ATAC
###################
###################

module load conda/anaconda3
source activate atacseq_2024

cd /Run

ls -1 *.fastq.gz > SE_files.txt

######### Trimming #########

lista=/Run/SE_files.txt

for file in `cat $lista`
do
        echo $file   
        trimmomatic SE -threads 5 -phred33 $file "`basename $file .fastq.gz`.trimmomatic_out.fastq.gz" ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done


date
echo "DONE"

# move trimmed fastq into /Run/Trimmed

#########


cd /Run/Trimmed

ls -1 *.trimmomatic_out.fastq.gz > trimmed_files.txt


######### Mapping #########

lista=/Run/Trimmed/trimmed_files.txt

for file in `cat $lista`
do
	echo $file 
    bwa mem -t 4 human_ref $file | samtools sort -@4 -o "`basename $file .trimmomatic_out.fastq.gz`.bam" - \
	&& samtools index "`basename $file .trimmomatic_out.fastq.gz`.bam" \
   	&& samtools flagstat "`basename $file .trimmomatic_out.fastq.gz`.bam" > "`basename $file .trimmomatic_out.fastq.gz`_map_stats.txt"
done

# move mapped fastq into /Run/Trimmed/Mapped

#########


cd /Run/Trimmed/Mapped

ls -1 *.bam > mapped_files.txt


######### remove Mito and strange Chr and mark duplicate #########

lista=/Run/Trimmed/Mapped/mapped_files.txt

for file in `cat $lista`
do
	echo $file 

	samtools view -bh $file chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > "`basename $file .bam`_SC_subset.bam"	\
   	 && picard MarkDuplicates I="`basename $file .bam`_SC_subset.bam" O="`basename $file .bam`_SC_subset_dedup.bam" M="`basename $file .bam`_markdup_metrics.txt"	\
   	 && samtools index "`basename $file .bam`_SC_subset_dedup.bam"	\
   	 && samtools flagstat "`basename $file .bam`_SC_subset_dedup.bam" > "`basename $file .bam`_SC_subset_dedup_map_stats.txt"	\
   	 && samtools idxstats  "`basename $file .bam`_SC_subset_dedup.bam"  >  "`basename $file .bam`.idxstats"

done

# move Filtered fastq into /Run/Trimmed/Mapped/Filtered

#########


cd /Run/Trimmed/Mapped/Filtered

ls -1 *_SC_subset_dedup.bam > filtered_files.txt


#########


######### remove duplicate e black regions #########

lista=/Run/Trimmed/Mapped/Filtered/filtered_files.txt

for file in `cat $lista`
do

	alignmentSieve -b $file -o "`basename $file .bam`_blackremove.bam" -p 10 --ignoreDuplicates --blackListFileName /reference/hg38-blacklist.v2.bed --minMappingQuality 10 --filterMetrics "`basename $file .bam`_shifted_log.txt"	\
	&& samtools sort "`basename $file .bam`_blackremove.bam" > "`basename $file .bam`_blackremove_sorted.bam"	\
	&& samtools index "`basename $file .bam`_blackremove_sorted.bam"

done

#########


######### in R: shift for Tn5 insertion #########

module load conda/anaconda3
source activate atacseq_2024

R
set.seed(123)


setwd("/Run/Trimmed/Mapped/Filtered")


library(ATACseqQC)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPpeakAnno)
library(Rsamtools)

# path to the processed bams file
Bam_1="1_DonorA_Gplus_Control_SC_subset_dedup_blackremove_sorted.bam"
Bam_2="2_DonorA_Gplus_NaCl_SC_subset_dedup_blackremove_sorted.bam"
Bam_3="3_DonorB_Gplus_Control_SC_subset_dedup_blackremove_sorted.bam"
Bam_4="4_DonorB_Gplus_NaCl_SC_subset_dedup_blackremove_sorted.bam"
Bam_5="5_DonorC_Gplus_Control_SC_subset_dedup_blackremove_sorted.bam"
Bam_6="6_DonorC_Gplus_NaCl_SC_subset_dedup_blackremove_sorted.bam"
Bam_7="7_DonorA_Gminus_Control_SC_subset_dedup_blackremove_sorted.bam"
Bam_8="8_DonorA_Gminus_NaCl_SC_subset_dedup_blackremove_sorted.bam"
Bam_9="9_DonorB_Gminus_Control_SC_subset_dedup_blackremove_sorted.bam"
Bam_10="10_DonorB_Gminus_NaCl_SC_subset_dedup_blackremove_sorted.bam"
Bam_11="11_DonorC_Gminus_Control_SC_subset_dedup_blackremove_sorted.bam"
Bam_12="12_DonorC_Gminus_NaCl_SC_subset_dedup_blackremove_sorted.bam"


## files will be saved into outPath respective to the working directory
outPath <- "Shifted"
dir.create(outPath)

possibleTag <- combn(LETTERS, 2)
possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]),
                 paste0(possibleTag[2, ], possibleTag[1, ]))



bamTop100_1 <- scanBam(BamFile(Bam_1, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags_1 <- names(bamTop100_1)[lengths(bamTop100_1)>0]

##

bamTop100_2 <- scanBam(BamFile(Bam_2, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags_2 <- names(bamTop100_2)[lengths(bamTop100_2)>0]

##

bamTop100_3 <- scanBam(BamFile(Bam_3, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags_3 <- names(bamTop100_3)[lengths(bamTop100_3)>0]

##

bamTop100_4 <- scanBam(BamFile(Bam_4, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags_4 <- names(bamTop100_4)[lengths(bamTop100_4)>0]

##

bamTop100_5 <- scanBam(BamFile(Bam_5, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags_5 <- names(bamTop100_5)[lengths(bamTop100_5)>0]

##

bamTop100_6 <- scanBam(BamFile(Bam_6, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags_6 <- names(bamTop100_6)[lengths(bamTop100_6)>0]

##

bamTop100_7 <- scanBam(BamFile(Bam_7, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags_7 <- names(bamTop100_7)[lengths(bamTop100_7)>0]

##

bamTop100_8 <- scanBam(BamFile(Bam_8, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags_8 <- names(bamTop100_8)[lengths(bamTop100_8)>0]

##

bamTop100_9 <- scanBam(BamFile(Bam_9, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags_9 <- names(bamTop100_9)[lengths(bamTop100_9)>0]

##

bamTop100_10 <- scanBam(BamFile(Bam_10, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags_10 <- names(bamTop100_10)[lengths(bamTop100_10)>0]

##

bamTop100_11 <- scanBam(BamFile(Bam_11, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags_11 <- names(bamTop100_11)[lengths(bamTop100_11)>0]

##

bamTop100_12 <- scanBam(BamFile(Bam_12, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags_12 <- names(bamTop100_12)[lengths(bamTop100_12)>0]

##


seqlev <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
which <- as(seqinfo(Hsapiens)[seqlev], "GRanges")

gal_1 <- readBamFile(Bam_1, tag=tags_1, which=which,asMates=FALSE, bigFile=TRUE)
gal_2 <- readBamFile(Bam_2, tag=tags_2, which=which,asMates=FALSE, bigFile=TRUE)
gal_3 <- readBamFile(Bam_3, tag=tags_3, which=which,asMates=FALSE, bigFile=TRUE)
gal_4 <- readBamFile(Bam_4, tag=tags_4, which=which,asMates=FALSE, bigFile=TRUE)
gal_5 <- readBamFile(Bam_5, tag=tags_5, which=which,asMates=FALSE, bigFile=TRUE)
gal_6 <- readBamFile(Bam_6, tag=tags_6, which=which,asMates=FALSE, bigFile=TRUE)
gal_7 <- readBamFile(Bam_7, tag=tags_7, which=which,asMates=FALSE, bigFile=TRUE)
gal_8 <- readBamFile(Bam_8, tag=tags_8, which=which,asMates=FALSE, bigFile=TRUE)
gal_9 <- readBamFile(Bam_9, tag=tags_9, which=which,asMates=FALSE, bigFile=TRUE)
gal_10 <- readBamFile(Bam_10, tag=tags_10, which=which,asMates=FALSE, bigFile=TRUE)
gal_11 <- readBamFile(Bam_11, tag=tags_11, which=which,asMates=FALSE, bigFile=TRUE)
gal_12 <- readBamFile(Bam_12, tag=tags_12, which=which,asMates=FALSE, bigFile=TRUE)


shiftedBamFile_1 <- file.path(outPath, "1_shifted.bam")
shiftedBamFile_2 <- file.path(outPath, "2_shifted.bam")
shiftedBamFile_3 <- file.path(outPath, "3_shifted.bam")
shiftedBamFile_4 <- file.path(outPath, "4_shifted.bam")
shiftedBamFile_5 <- file.path(outPath, "5_shifted.bam")
shiftedBamFile_6 <- file.path(outPath, "6_shifted.bam")
shiftedBamFile_7 <- file.path(outPath, "7_shifted.bam")
shiftedBamFile_8 <- file.path(outPath, "8_shifted.bam")
shiftedBamFile_9 <- file.path(outPath, "9_shifted.bam")
shiftedBamFile_10 <- file.path(outPath, "10_shifted.bam")
shiftedBamFile_11 <- file.path(outPath, "11_shifted.bam")
shiftedBamFile_12 <- file.path(outPath, "12_shifted.bam")


galc <- shiftGAlignments(gal_1, positive = 4L, negative = 5L, outbam=shiftedBamFile_1)
galm <- shiftGAlignments(gal_2, positive = 4L, negative = 5L, outbam=shiftedBamFile_2)
galn <- shiftGAlignments(gal_3, positive = 4L, negative = 5L, outbam=shiftedBamFile_3)
galo <- shiftGAlignments(gal_4, positive = 4L, negative = 5L, outbam=shiftedBamFile_4)
galp <- shiftGAlignments(gal_5, positive = 4L, negative = 5L, outbam=shiftedBamFile_5)
galq <- shiftGAlignments(gal_6, positive = 4L, negative = 5L, outbam=shiftedBamFile_6)
galr <- shiftGAlignments(gal_7, positive = 4L, negative = 5L, outbam=shiftedBamFile_7)
gals <- shiftGAlignments(gal_8, positive = 4L, negative = 5L, outbam=shiftedBamFile_8)
galt <- shiftGAlignments(gal_9, positive = 4L, negative = 5L, outbam=shiftedBamFile_9)
gala <- shiftGAlignments(gal_10,  positive = 4L, negative = 5L, outbam=shiftedBamFile_10)
galb <- shiftGAlignments(gal_11, positive = 4L, negative = 5L, outbam=shiftedBamFile_11)
gald <- shiftGAlignments(gal_12, positive = 4L, negative = 5L, outbam=shiftedBamFile_12)


#########


cd /Run/Trimmed/Mapped/Filtered/Shifted
ls -1 *.bam > shifted_files.txt


#########


######### call peaks MACS3 #########

lista=/Run/Trimmed/Mapped/Filtered/Shifted/shifted_files.txt

for file in `cat $lista`
do
	echo $file 
	macs3 callpeak -t $file -n "`basename $file .bam`_" -f BAM --outdir Run/Peaks/ -g hs --nomodel --shift -100 --extsize 200 --keep-dup all -B
done

#########


cd /Run/Peaks
ls -1 *__treat_pileup.bdg > track_files.txt


#########

# FIGURE 5K
######### Creating a browser track so we can look at the peaks in the UCSC Genome Browser (FIG5_K) #########

lista=/Run/Peaks/track_files.txt

for file in `cat $lista`
do
	LC_COLLATE=C sort -k1,1 -k2,2n $file > "`basename $file .bdg`_sorted.bdg"	\
	&& bedGraphToBigWig "`basename $file .bdg`_sorted.bdg" /reference/hg38.chrom.sizes "`basename $file __treat_pileup.bdg`_peaks.bw"
done

#########


cd /Run/Peaks
ls -1 *__peaks.narrowPeak > peaks_files.txt


#########


######### Remove blacklist regions from Peaks #########


lista=/Run/Peaks/peaks_files.txt

for file in `cat $lista`
do

	bedtools intersect -v -a $file  -b /reference/hg38-blacklist.v2.bed  >| "`basename $file __peaks.narrowPeak`_BlackFiltered_peaks.narrowPeak"

done


#########


######### Consensus peaks in at least 2 samples per condition #########


# Concatenate narrowPeak files, coordinate sort, then merge peaks within 10 bp
cd /Run/Peaks/consensus_at_least2

## Gplus_Control ## 

cat /Run/Peaks/1_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/3_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/5_shifted_BlackFiltered_peaks.narrowPeak > /Run/Peaks/consensus_at_least2/Gplus_Control_peaks_123
sort -k1,1 -k2,2n -k3,3n /Run/Peaks/consensus_at_least2/Gplus_Control_peaks_123 > /Run/Peaks/consensus_at_least2/Gplus_Control_peaks_123_sort
cat /Run/Peaks/consensus_at_least2/Gplus_Control_peaks_123_sort | awk '{print $1"\t"$2"\t"$3"\t"$7}' > /Run/Peaks/consensus_at_least2/Gplus_Control_peaks_123_sort.bed
bedtools merge -d 10 -c 4 -o mean -i /Run/Peaks/consensus_at_least2/Gplus_Control_peaks_123_sort.bed > /Run/Peaks/consensus_at_least2/Gplus_Control_peaks_mergedpeaks.bed

wc -l Gplus_Control_peaks_mergedpeaks.bed
77483 Gplus_Control_peaks_mergedpeaks.bed

# Write out bed file with showing the number of replicates that support peaks in summary file

bedtools intersect -wa -c -a Gplus_Control_peaks_mergedpeaks.bed -b /Run/Peaks/1_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/3_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/5_shifted_BlackFiltered_peaks.narrowPeak -sorted -F 1.0 > Gplus_Control_mergedpeak_replicates.bed

# Filter file for peaks that have greater than or equal to 2 samples, sort, then remove the column with the number of replicates:
awk '$5 >= 2 {print}' Gplus_Control_mergedpeak_replicates.bed > Gplus_Control_mergedpeak_replicates_filter.bed
sort -k1,1 -k2,2n -k3,3n Gplus_Control_mergedpeak_replicates_filter.bed > Gplus_Control_mergedpeak_replicates_filter_sort.bed
cat Gplus_Control_mergedpeak_replicates_filter_sort.bed | awk '{print $1"\t"$2"\t"$3}' > Gplus_Control_peak_replicates_atleast2.bed

wc -l Gplus_Control_peak_replicates_atleast2.bed
56801 Gplus_Control_peak_replicates_atleast2.bed


#########


## Gplus_NaCl ## 

cat /Run/Peaks/2_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/4_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/6_shifted_BlackFiltered_peaks.narrowPeak > /Run/Peaks/consensus_at_least2/Gplus_NaCl_peaks_123
sort -k1,1 -k2,2n -k3,3n /Run/Peaks/consensus_at_least2/Gplus_NaCl_peaks_123 > /Run/Peaks/consensus_at_least2/Gplus_NaCl_peaks_123_sort
cat /Run/Peaks/consensus_at_least2/Gplus_NaCl_peaks_123_sort | awk '{print $1"\t"$2"\t"$3"\t"$7}' > /Run/Peaks/consensus_at_least2/Gplus_NaCl_peaks_123_sort.bed
bedtools merge -d 10 -c 4 -o mean -i /Run/Peaks/consensus_at_least2/Gplus_NaCl_peaks_123_sort.bed > /Run/Peaks/consensus_at_least2/Gplus_NaCl_peaks_mergedpeaks.bed

wc -l Gplus_NaCl_peaks_mergedpeaks.bed
135488 Gplus_NaCl_peaks_mergedpeaks.bed

# Write out bed file with showing the number of replicates that support peaks in summary file

bedtools intersect -wa -c -a Gplus_NaCl_peaks_mergedpeaks.bed -b /Run/Peaks/2_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/4_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/6_shifted_BlackFiltered_peaks.narrowPeak -sorted -F 1.0 > Gplus_NaCl_mergedpeak_replicates.bed

# Filter file for peaks that have greater than or equal to 2 samples, sort, then remove the column with the number of replicates:
awk '$5 >= 2 {print}' Gplus_NaCl_mergedpeak_replicates.bed > Gplus_NaCl_mergedpeak_replicates_filter.bed
sort -k1,1 -k2,2n -k3,3n Gplus_NaCl_mergedpeak_replicates_filter.bed > Gplus_NaCl_mergedpeak_replicates_filter_sort.bed
cat Gplus_NaCl_mergedpeak_replicates_filter_sort.bed | awk '{print $1"\t"$2"\t"$3}' > Gplus_NaCl_peak_replicates_atleast2.bed

wc -l Gplus_NaCl_peak_replicates_atleast2.bed
100765 Gplus_NaCl_peak_replicates_atleast2.bed


#########


## Gminus_Control ## 

cat /Run/Peaks/7_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/9_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/11_shifted_BlackFiltered_peaks.narrowPeak > /Run/Peaks/consensus_at_least2/Gminus_Control_peaks_123
sort -k1,1 -k2,2n -k3,3n /Run/Peaks/consensus_at_least2/Gminus_Control_peaks_123 > /Run/Peaks/consensus_at_least2/Gminus_Control_peaks_123_sort
cat /Run/Peaks/consensus_at_least2/Gminus_Control_peaks_123_sort | awk '{print $1"\t"$2"\t"$3"\t"$7}' > /Run/Peaks/consensus_at_least2/Gminus_Control_peaks_123_sort.bed
bedtools merge -d 10 -c 4 -o mean -i /Run/Peaks/consensus_at_least2/Gminus_Control_peaks_123_sort.bed > /Run/Peaks/consensus_at_least2/Gminus_Control_peaks_mergedpeaks.bed

wc -l Gminus_Control_peaks_mergedpeaks.bed
86274 Gminus_Control_peaks_mergedpeaks.bed

# Write out bed file with showing the number of replicates that support peaks in summary file

bedtools intersect -wa -c -a Gminus_Control_peaks_mergedpeaks.bed -b /Run/Peaks/7_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/9_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/11_shifted_BlackFiltered_peaks.narrowPeak -sorted -F 1.0 > Gminus_Control_mergedpeak_replicates.bed

# Filter file for peaks that have greater than or equal to 2 samples, sort, then remove the column with the number of replicates:
awk '$5 >= 2 {print}' Gminus_Control_mergedpeak_replicates.bed > Gminus_Control_mergedpeak_replicates_filter.bed
sort -k1,1 -k2,2n -k3,3n Gminus_Control_mergedpeak_replicates_filter.bed > Gminus_Control_mergedpeak_replicates_filter_sort.bed
cat Gminus_Control_mergedpeak_replicates_filter_sort.bed | awk '{print $1"\t"$2"\t"$3}' > Gminus_Control_peak_replicates_atleast2.bed

wc -l Gminus_Control_peak_replicates_atleast2.bed
67959 Gminus_Control_peak_replicates_atleast2.bed


#########


## Gminus_NaCl ## 

cat /Run/Peaks/8_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/10_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/12_shifted_BlackFiltered_peaks.narrowPeak > /Run/Peaks/consensus_at_least2/Gminus_NaCl_peaks_123
sort -k1,1 -k2,2n -k3,3n /Run/Peaks/consensus_at_least2/Gminus_NaCl_peaks_123 > /Run/Peaks/consensus_at_least2/Gminus_NaCl_peaks_123_sort
cat /Run/Peaks/consensus_at_least2/Gminus_NaCl_peaks_123_sort | awk '{print $1"\t"$2"\t"$3"\t"$7}' > /Run/Peaks/consensus_at_least2/Gminus_NaCl_peaks_123_sort.bed
bedtools merge -d 10 -c 4 -o mean -i /Run/Peaks/consensus_at_least2/Gminus_NaCl_peaks_123_sort.bed > /Run/Peaks/consensus_at_least2/Gminus_NaCl_peaks_mergedpeaks.bed

wc -l Gminus_NaCl_peaks_mergedpeaks.bed
136008 Gminus_NaCl_peaks_mergedpeaks.bed

# Write out bed file with showing the number of replicates that support peaks in summary file

bedtools intersect -wa -c -a Gminus_NaCl_peaks_mergedpeaks.bed -b /Run/Peaks/8_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/10_shifted_BlackFiltered_peaks.narrowPeak /Run/Peaks/12_shifted_BlackFiltered_peaks.narrowPeak -sorted -F 1.0 > Gminus_NaCl_mergedpeak_replicates.bed

# Filter file for peaks that have greater than or equal to 2 samples, sort, then remove the column with the number of replicates:
awk '$5 >= 2 {print}' Gminus_NaCl_mergedpeak_replicates.bed > Gminus_NaCl_mergedpeak_replicates_filter.bed
sort -k1,1 -k2,2n -k3,3n Gminus_NaCl_mergedpeak_replicates_filter.bed > Gminus_NaCl_mergedpeak_replicates_filter_sort.bed
cat Gminus_NaCl_mergedpeak_replicates_filter_sort.bed | awk '{print $1"\t"$2"\t"$3}' > Gminus_NaCl_peak_replicates_atleast2.bed

wc -l Gminus_NaCl_peak_replicates_atleast2.bed
101131 Gminus_NaCl_peak_replicates_atleast2.bed


#########


######### in R: Identifying a set of non-redundant Peaks #########


module load conda/anaconda3
source activate atacseq_2024

R
set.seed(123)


setwd("/Run/Peaks/Analisys")

peaks <- dir("/Run/Peaks/consensus_at_least2/", pattern = "*replicates_atleast2.bed", full.names = TRUE)
myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)

library(purrr)
library(GenomicRanges)

allPeaksSet_nR <- reduce(unlist(GRangesList(myPeaks)))
overlap <- list()
for (i in 1:length(myPeaks)) {
    overlap[[i]] <- allPeaksSet_nR %over% myPeaks[[i]]
}
overlapMatrix <- do.call(cbind, overlap)
colnames(overlapMatrix) <- basename(peaks)
mcols(allPeaksSet_nR) <- overlapMatrix

x <- paste0("all_condition_merged_macs3_", 1:length(allPeaksSet_nR))

allPeaksSet_nR <- setNames(allPeaksSet_nR, x)
length(allPeaksSet_nR)

library(rtracklayer)

to_bed <- data.frame(seqnames=seqnames(allPeaksSet_nR), starts=start(allPeaksSet_nR)-1, ends=end(allPeaksSet_nR))
write.table(to_bed, file="all_condition_merged_peaks_least2.bed", quote=F, sep="\t", row.names=F, col.names=F)

to_saf <- data.frame(names(allPeaksSet_nR),seqnames=seqnames(allPeaksSet_nR), starts=start(allPeaksSet_nR)-1, ends=end(allPeaksSet_nR),strands=c(rep(".", length(allPeaksSet_nR))))
write.table(to_saf, file="all_condition_merged_peaks_least2.saf", quote=F, sep="\t", row.names=F, col.names=F)

to_id <- data.frame(seqnames=seqnames(allPeaksSet_nR), starts=start(allPeaksSet_nR)-1, ends=end(allPeaksSet_nR),names(allPeaksSet_nR),scores=c(rep("0", length(allPeaksSet_nR))),strands=c(rep(".", length(allPeaksSet_nR))))
write.table(to_id, file="all_condition_merged_peaks_least2_peaksID.bed", quote=F, sep="\t", row.names=F, col.names=F)


#########


######### Obtain the Raw matrix #########


cd /Run/Peaks/Analisys

featureCounts -T 10 -F SAF -a /Run/Peaks/Analisys/all_condition_merged_peaks_least2.saf --fracOverlap 0.2 -o /Run/Peaks/Analisys/all_condition_merged_peaks_atleast2.counts /Run/Trimmed/Mapped/Filtered/Shifted/1_shifted.bam /Run/Trimmed/Mapped/Filtered/Shifted/3_shifted.bam /Run/Trimmed/Mapped/Filtered/Shifted/5_shifted.bam /Run/Trimmed/Mapped/Filtered/Shifted/2_shifted.bam /Run/Trimmed/Mapped/Filtered/Shifted/4_shifted.bam /Run/Trimmed/Mapped/Filtered/Shifted/6_shifted.bam /Run/Trimmed/Mapped/Filtered/Shifted/7_shifted.bam /Run/Trimmed/Mapped/Filtered/Shifted/9_shifted.bam /Run/Trimmed/Mapped/Filtered/Shifted/11_shifted.bam /Run/Trimmed/Mapped/Filtered/Shifted/8_shifted.bam /Run/Trimmed/Mapped/Filtered/Shifted/10_shifted.bam /Run/Trimmed/Mapped/Filtered/Shifted/12_shifted.bam


# We should remove the first line starting with #, as it can interfere with the way R reads in data:
cd /Run/Peaks/Analisys
awk '(NR>1)' all_condition_merged_peaks_atleast2.counts > all_condition_merged_peaks_atleast2.counts.tsv


#########


######### in R: correlation #########


module load conda/anaconda3
source activate atacseq_2024

R
set.seed(123)


setwd("/Run/Peaks/Analisys")


library(DiffBind)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

### Load the peak files from MACS to DiffBind

# Gplus_Control
experiment <- dba.peakset(NULL,
                      	peaks="/Run/Peaks/Analisys/all_condition_merged_peaks_least2.bed",
                      	peak.caller="bed", sampID="DonorA_Gplus_Control",condition = "Gplus_Control", replicate = 1, bamReads = "/Run/Trimmed/Mapped/Filtered/Shifted/1_shifted.bam")

experiment <- dba.peakset(experiment,
                      	peaks="/Run/Peaks/Analisys/all_condition_merged_peaks_least2.bed",
                      	peak.caller="bed", sampID="DonorB_Gplus_Control",condition = "Gplus_Control", replicate = 2, bamReads = "/Run/Trimmed/Mapped/Filtered/Shifted/3_shifted.bam")

experiment <- dba.peakset(experiment,
                      	peaks="/Run/Peaks/Analisys/all_condition_merged_peaks_least2.bed",
                      	peak.caller="bed", sampID="DonorC_Gplus_Control",condition = "Gplus_Control", replicate = 3, bamReads = "/Run/Trimmed/Mapped/Filtered/Shifted/5_shifted.bam")


# Gplus_NaCl
experiment <- dba.peakset(experiment,
                      	peaks="/Run/Peaks/Analisys/all_condition_merged_peaks_least2.bed",
                      	peak.caller="bed", sampID="DonorA_Gplus_NaCl",condition = "Gplus_NaCl", replicate = 1, bamReads = "/Run/Trimmed/Mapped/Filtered/Shifted/2_shifted.bam")

experiment <- dba.peakset(experiment,
                      	peaks="/Run/Peaks/Analisys/all_condition_merged_peaks_least2.bed",
                      	peak.caller="bed", sampID="DonorB_Gplus_NaCl",condition = "Gplus_NaCl", replicate = 2, bamReads = "/Run/Trimmed/Mapped/Filtered/Shifted/4_shifted.bam")

experiment <- dba.peakset(experiment,
                      	peaks="/Run/Peaks/Analisys/all_condition_merged_peaks_least2.bed",
                      	peak.caller="bed", sampID="DonorC_Gplus_NaCl",condition = "Gplus_NaCl", replicate = 3, bamReads = "/Run/Trimmed/Mapped/Filtered/Shifted/6_shifted.bam")


# Gminus_Control
experiment <- dba.peakset(experiment,
                      	peaks="/Run/Peaks/Analisys/all_condition_merged_peaks_least2.bed",
                      	peak.caller="bed", sampID="DonorA_Gminus_Control",condition = "Gminus_Control", replicate = 1, bamReads = "/Run/Trimmed/Mapped/Filtered/Shifted/7_shifted.bam")

experiment <- dba.peakset(experiment,
                      	peaks="/Run/Peaks/Analisys/all_condition_merged_peaks_least2.bed",
                      	peak.caller="bed", sampID="DonorB_Gminus_Control",condition = "Gminus_Control", replicate = 2, bamReads = "/Run/Trimmed/Mapped/Filtered/Shifted/9_shifted.bam")

experiment <- dba.peakset(experiment,
                      	peaks="/Run/Peaks/Analisys/all_condition_merged_peaks_least2.bed",
                      	peak.caller="bed", sampID="DonorC_Gminus_Control",condition = "Gminus_Control", replicate = 3, bamReads = "/Run/Trimmed/Mapped/Filtered/Shifted/11_shifted.bam")


# Gminus_NaCl
experiment <- dba.peakset(experiment,
                      	peaks="/Run/Peaks/Analisys/all_condition_merged_peaks_least2.bed",
                      	peak.caller="bed", sampID="DonorA_Gminus_NaCl",condition = "Gminus_NaCl", replicate = 1, bamReads = "/Run/Trimmed/Mapped/Filtered/Shifted/8_shifted.bam")

experiment <- dba.peakset(experiment,
                      	peaks="/Run/Peaks/Analisys/all_condition_merged_peaks_least2.bed",
                      	peak.caller="bed", sampID="DonorB_Gminus_NaCl",condition = "Gminus_NaCl", replicate = 2, bamReads = "/Run/Trimmed/Mapped/Filtered/Shifted/10_shifted.bam")

experiment <- dba.peakset(experiment,
                      	peaks="/Run/Peaks/Analisys/all_condition_merged_peaks_least2.bed",
                      	peak.caller="bed", sampID="DonorC_Gminus_NaCl",condition = "Gminus_NaCl", replicate = 3, bamReads = "/Run/Trimmed/Mapped/Filtered/Shifted/12_shifted.bam")



#### 
experiment_counts <- dba.count(experiment, bParallel = TRUE, score = DBA_SCORE_READS, bRemoveDuplicates = FALSE)

experiment_counts$config$RunParallel <- TRUE
profiles <- dba.plotProfile(experiment_counts, maxSites = 50000, distanceAround=1000)

# FIGURE 5J

pdf("FIG5_J.pdf", width = 30, height = 30)
plot(experiment_counts, main="Correlation plot of D5 samples")
dev.off()

data <- plot(experiment_counts, main="Correlation plot of D5 samples")
write.table(data, "raw_data_FIG5_J.txt", sep="\t", col.name = TRUE, quote = FALSE)


#########


######### in R: Differential analysis Gplus_Control VS Gplus_NaCl #########


module load conda/anaconda3
source activate atacseq_2024

R
set.seed(123)

setwd("/Run/Peaks/Analisys")

library(edgeR)
library(EDASeq)
library(GenomicAlignments)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(wesanderson)
library(Hmisc)
library(dplyr)
library(ggplot2)
library(ChIPseeker)
library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# After peak calling one may want to visualise distribution of peaks locations over the whole genome. Function covplot calculates coverage of peaks regions over chromosomes.
pth2peaks_bed <- "/Run/Peaks/Analisys/all_condition_merged_peaks_least2_peaksID.bed"

peaks.bed <- read.table(pth2peaks_bed, sep="\t", header=FALSE, blank.lines.skip=TRUE)
rownames(peaks.bed) <- peaks.bed[,4]

peaks.gr <- GRanges(seqnames=peaks.bed[,1], ranges=IRanges(peaks.bed[,2], peaks.bed[,3]), strand="*", mcols=data.frame(peakID=peaks.bed[,4]))

# annotate peaks with closest genomic features
bed.annot <- annotatePeak(peaks.gr, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
bed.annot

annot_peaks <- as.data.frame(bed.annot)

ff <- FaFile("/reference/WholeGenomeFasta/genome.fa")

# We can read in the data, format it and define experimental groups:
cnt_table <- read.table("/Run/Peaks/Analisys/all_condition_merged_peaks_atleast2.counts.tsv", sep="\t", header=TRUE, blank.lines.skip=TRUE)

rownames(cnt_table) <- cnt_table$Geneid

# update colnames of this count table
colnames(cnt_table) <- c("Geneid","Chr","Start","End","Strand","Length","DonorA_Gplus_Control","DonorB_Gplus_Control","DonorC_Gplus_Control","DonorA_Gplus_NaCl","DonorB_Gplus_NaCl","DonorC_Gplus_NaCl","DonorA_Gminus_Control","DonorB_Gminus_Control","DonorC_Gminus_Control","DonorA_Gminus_NaCl","DonorB_Gminus_NaCl","DonorC_Gminus_NaCl")

groups <- factor(c(rep("Gplus_Control",3),rep("Gplus_NaCl",3),rep("Gminus_Control",3),rep("Gminus_NaCl",3)))

# this data frame contains only read counts to peaks on assembled chromosomes
reads.peak <- cnt_table[,c(7:18)]

reads.peak <- as.matrix(reads.peak)

# DGE object
dge_data <- DGEList(reads.peak, group=groups)

keep <- filterByExpr(dge_data, group=groups)

dge_data <- dge_data[keep,,keep.lib.sizes=FALSE]


dge_data <- calcNormFactors(dge_data)
dge_data$samples

esperimento <- as.data.frame(cbind(colnames(dge_data),c(rep(c("DonorA","DonorB","DonorC"),4)),c(rep("Gplus_Control",3),rep("Gplus_NaCl",3),rep("Gminus_Control",3),rep("Gminus_NaCl",3))))
colnames(esperimento) <- c("Sample_ID", "Donor_ID", "Treatment")
esperimento


col_use <- c(rep("green",3),rep("black",3),rep("red",3),rep("blue",3))

design <- model.matrix(~0+groups)
dge_data <- estimateDisp(dge_data,design,robust=TRUE)

# calcola i CPM
matrice_cpm <- cpm(dge_data)

matrice_cpm_anno <- merge(annot_peaks,matrice_cpm,by.x="mcols.peakID",by.y="row.names",sort=FALSE)

#calcola i log2CPM
matrice_log2cpm <- log2(matrice_cpm+1)

matrice_log2cpm_anno <- merge(annot_peaks,matrice_log2cpm,by.x="mcols.peakID",by.y="row.names",sort=FALSE)

library(gplots)

fit <- glmQLFit(dge_data,design)


con_Gplus_Control_vs_Gplus_NaCl <- makeContrasts(groupsGplus_NaCl - groupsGplus_Control, levels=design)
qlf_contVStrat <- glmQLFTest(fit, contrast = con_Gplus_Control_vs_Gplus_NaCl)

summary(decideTests(qlf_contVStrat))

topTags(qlf_contVStrat)
GeniDE_qlf_contVStrat <- topTags(qlf_contVStrat, n = dim(qlf_contVStrat)[1])
toptag <- GeniDE_qlf_contVStrat@.Data[[1]]

DEG_FDR005 <- toptag[toptag$FDR<0.05,]

matrice_log2cpm_anno_new <- matrice_log2cpm_anno[,-c(25:30)]


# Let’s add more peak information
DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno <- merge(DEG_FDR005,matrice_log2cpm_anno_new,by.x="row.names",by.y="mcols.peakID",sort=FALSE)


# filter fdr<0.05 and Promoter (<=1kb)
DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno_prom <- DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno[DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno$annotation %in% "Promoter (<=1kb)",]

# downreg
nrow(DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno_prom[DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno_prom$logFC<0,])

# upreg
nrow(DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno_prom[DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno_prom$logFC>0,])

DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno_prom$unique_symbol <- paste(DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno_prom$SYMBOL,rownames(DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno_prom), sep = "_")



# Make Homer input
# filter only promoter (+- 1000 TSS) and FDR < 0.05

########
# downreg
prom_fc_neg <- DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno_prom[DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno_prom$logFC<0,]

prom_peaks_Gplus_Control_homer <- prom_fc_neg[,c("seqnames","start","end","Row.names")]
prom_peaks_Gplus_Control_homer$v1 <- "0"
prom_peaks_Gplus_Control_homer$v2 <- "+"
prom_peaks_Gplus_Control_homer <- prom_peaks_Gplus_Control_homer[order(prom_peaks_Gplus_Control_homer$Row.names),]
write.table(prom_peaks_Gplus_Control_homer, file="Gplus_Control_vs_Gplus_NaCl_promoter_peaks_homer_Gplus_Control.bed", quote=F, sep="\t", row.names=FALSE, col.names=FALSE)


########
# upreg
prom_fc_pos <- DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno_prom[DEG_FDR005_Gplus_Control_vs_Gplus_NaCl_anno_prom$logFC>0,]

prom_peaks_Gplus_NaCl_homer <- prom_fc_pos[,c("seqnames","start","end","Row.names")]
prom_peaks_Gplus_NaCl_homer$v1 <- "0"
prom_peaks_Gplus_NaCl_homer$v2 <- "+"
prom_peaks_Gplus_NaCl_homer <- prom_peaks_Gplus_NaCl_homer[order(prom_peaks_Gplus_NaCl_homer$Row.names),]
write.table(prom_peaks_Gplus_NaCl_homer, file="Gplus_Control_vs_Gplus_NaCl_promoter_peaks_homer_Gplus_NaCl.bed", quote=F, sep="\t", row.names=FALSE, col.names=FALSE)


#########


######### Convert from Bed to Fasta #########


cd /Run/Peaks/Analisys

bedtools getfasta -fi /reference/WholeGenomeFasta/genome.fa -bed /Run/Peaks/Analisys/Gplus_Control_vs_Gplus_NaCl_promoter_peaks_homer_Gplus_Control.bed -fo /Run/Peaks/Analisys/Gplus_Control_vs_Gplus_NaCl_promoter_peaks_homer_Gplus_Control.fa

bedtools getfasta -fi /reference/WholeGenomeFasta/genome.fa -bed /Run/Peaks/Analisys/Gplus_Control_vs_Gplus_NaCl_promoter_peaks_homer_Gplus_NaCl.bed -fo /Run/Peaks/Analisys/Gplus_Control_vs_Gplus_NaCl_promoter_peaks_homer_Gplus_NaCl.fa


#########


######### Homer TFs Analysis (FIG5_I) #########

# FIGURE 5I

findMotifs.pl /Run/Peaks/Analisys/Gplus_Control_vs_Gplus_NaCl_promoter_peaks_homer_Gplus_NaCl.fa fasta /Run/Peaks/Analisys/Homer/DAR_Gplus_NaCl_promoter_peaks_WithBackground -fasta /Run/Peaks/Analisys/Gplus_Control_vs_Gplus_NaCl_promoter_peaks_homer_Gplus_Control.fa -p 10

findMotifs.pl /Run/Peaks/Analisys/Gplus_Control_vs_Gplus_NaCl_promoter_peaks_homer_Gplus_Control.fa fasta /Run/Peaks/Analisys/Homer/DAR_Gplus_Control_promoter_peaks_WithBackground -fasta /Run/Peaks/Analisys/Gplus_Control_vs_Gplus_NaCl_promoter_peaks_homer_Gplus_NaCl.fa -p 10


#########
