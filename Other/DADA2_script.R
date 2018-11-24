library(dada2); packageVersion("dada2")

path="C:/My Directory"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs = sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs = sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])

plotQualityProfile(fnRs[1:2])

filt_path = file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs = file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs = file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

errF = learnErrors(filtFs, multithread=TRUE)
errR = learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

derepFs = derepFastq(filtFs, verbose=TRUE)
derepRs = derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) = sample.names
names(derepRs) = sample.names

dadaFs = dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab = makeSequenceTable(mergers)

dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
