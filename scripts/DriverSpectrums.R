ParseAnnovarTranscript <- function(string){
  annotations <- strsplit(string, ',')[[1]]
  gene <- paste(sapply(annotations, function(x) strsplit(x, ':')[[1]][2]), collapse = ':')
  return(gene)
}

ParseAnnovarMutation <- function(string){
  annotations <- strsplit(string, ',')[[1]]
  mutation <- paste(sapply(annotations, function(x) substr(strsplit(x, ':')[[1]][5], 3, nchar(strsplit(x, ':')[[1]][4]))), collapse = ':')
}

ParseAnnovarGene <- function(string){
  mutations <- strsplit(string, ',')[[1]]
  gene <- paste(unique(sapply(mutations, function(x) strsplit(x, ':')[[1]][1])), collapse = ':')
  return(gene)
}

####################################################################################################################################


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(annotate)
library(org.Hs.eg.db)
require('SomaticSignatures')
# load mutation data

# read in oncogenes and assign ensemble gene.id column
setwd("/Users/[user]/data/ref/")
genes <- read.table('Vogelstein.Cancer.Genome.Landscapes.Table.S2a.edit.2.Aug.2017.csv', sep = ',', header = T, stringsAsFactors = F)
map <- read.table('ensembl.id.to.hgnc.symbol.txt', sep = '\t', header = T, stringsAsFactors = F)
genes$name <- genes$gene_name
genes$name[match(c('FAM123B', 'MLL2', 'MLL3'), genes$gene_name)] <- c('AMER1', 'KMT2D', 'KMT2C')
# assigne ensembl.ids to ogs
genes$ensembl.id <- NA
for (i in 1:nrow(genes)){
  idx <- match(genes$name[i], map$HGNC.symbol)
  if (!is.na(idx)){
    genes$ensembl.id[i] <- map$Ensembl.Gene.ID[idx]
  }
}
ogs <- genes[which(genes$classification == 'Oncogene'), ]
eg.map <- as.list(org.Hs.egSYMBOL2EG)
genes$entrez.id <- as.character(eg.map[genes$gene_name])

canonical <- read.table('/Users/[user]/data/ref/knownCanonical.5.June.2017.txt', col.names = c('chr', 'start', 'end', 'row', 'transcript.id', 'Ensembl.Gene.ID.full'))
canonical$Ensembl.Gene.ID <- sapply(canonical$Ensembl.Gene.ID.full, function(x) strsplit(as.character(x), '\\.')[[1]][1])
genes$transcript.id <- canonical$transcript.id[match(genes$ensembl.id, canonical$Ensembl.Gene.ID)]
cross.ref <- read.table('/Users/[user]/data/ref/knownToRefSeq.6.June.tsv', col.names = c('transcript.id', 'transcript'))
genes$transcript <- cross.ref$transcript[match(genes$transcript.id, cross.ref$transcript.id)]
genes <- genes[which(!(genes$name == 'U2AF1')),]

# read in metadata
msl <- read.table("/Users/[user]/data/driver_spectrum/metadata/MasterSampleList.csv", sep = ',', header = T, stringsAsFactors = F)
cancer.groupings <- read.table("/Users/[user]/data/driver_spectrum/metadata/Disease.Classification.csv", sep = ',', header = T, stringsAsFactors = F)
cancer.groupings$disease.group <- sapply(1:nrow(cancer.groupings), function(x) if(!is.na(cancer.groupings$Signature.Group[x])) {cancer.groupings$Signature.Group[x]} else {cancer.groupings$Disease[x]})
tcga.samples <- read.table("/Users/[user]/data/driver_spectrum/metadata/TCGA.Files.csv", sep = ',', header = T, stringsAsFactors = F)
tcga.samples$Data.ID <- sapply(tcga.samples$Paths, function(x) strsplit(x, '/')[[1]][1])
# read in data
short.names <- msl$Study.ID

# read in Alexandrov singatures and normalise to exome
prc <- read.table('/Users/[user]/data/ref/signatures_probabilities.txt', sep = '\t', header = T, stringsAsFactors = F)
unordered.processes <- paste(gsub('>', '', prc$Substitution.Type), paste(substr(prc$Trinucleotide, 1, 1), substr(prc$Trinucleotide, 3, 3), sep = '.'))
prc <- prc[order(unordered.processes), ]
processes <- sort(unordered.processes)
tnf <- read.table("/Users/[user]/data/ref/trinucleotide.frequencies.tsv", header = T, stringsAsFactors = F)
cnf <- data.frame(process = processes, type = as.character(sapply(processes, function(x) paste(substr(x,4,4), substr(x,1,1), substr(x,6,6), sep = ''))))
cnf$genome <- tnf$genome[match(cnf$type, tnf$type)]
cnf$exome <- tnf$exome[match(cnf$type, tnf$type)]
prc.proven <- prc[, paste('Signature.', 1:30, sep = '')]
# normalise signatures by trinucleotide frequencies in the exome
prc.norm <- prc
for(i in grep('Signature', colnames(prc))){
  prc.norm[, i] <- prc[, i] * (cnf$exome / cnf$genome)
  prc.norm[, i] <- prc.norm[, i] / sum(prc.norm[, i])
}
prc.proven.exome <- prc.norm[, paste('Signature.', 1:30, sep = '')]
# normalise signatures by equal trinculeotide frequencies
prc.norm <- prc
for(i in grep('Signature', colnames(prc))){
  prc.norm[, i] <- prc[, i] * (1 / cnf$genome)
  prc.norm[, i] <- prc.norm[, i] / sum(prc.norm[, i])
}
prc.proven.even <- prc.norm[, paste('Signature.', 1:30, sep = '')]

genome.signature.scaling.factors <- apply(prc.proven.even, 2, function(x) sum(x * cnf$genome))
exome.signature.scaling.factors <- apply(prc.proven.even, 2, function(x) sum(x * cnf$exome))

####################################################################################################################################
options(warn = 2)

# read in data from all files [VERY SLOW]
icgc.errors <- c()
error.names <- c()
for (short.name in sub){
  if (msl$ICGC.portal[match(short.name, msl$Study.ID)] == 'no'){
    next
  }
  file <- paste("/Users/[user]/apocrita/data_mount_point_2/daniel_import/signatures/ICGC/ssm.", msl$ICGC.Study.Abbr[match(short.name, msl$Study.ID)], ".tsv", sep = '')
  
  print(paste('Parsing', file))
  # read in data
  data <- read.table(file = file, sep = '\t', header = T, stringsAsFactors = F)
  
  # SNVs only
  data <- data[which(data$mutation_type == 'single base substitution'), ]
  
  # remove cases where mutated_from_allele equals mutated_to_allel
  if (length(which(data$mutated_from_allele == data$mutated_to_allele)) > 0){
    data <- data[-which(data$mutated_from_allele == data$mutated_to_allele), ]
  }
  
  # remove duplicates
  data$mutation <- paste(data$icgc_donor_id, data$chromosome, data$chromosome_start)
  data <- data[unique(match(data$mutation, data$mutation)),]
  
  # remove mitochondrial mutations
  data <- data[which(data$chromosome %in% c(1:22, 'X', 'Y')), ]
  
  # save annotations
  annotation <- data.frame(chromosome = data$chromosome,  pos = data$chromosome_start, 
                           ref = data$mutated_from_allele, alt = data$mutated_to_allele,
                           gene = data$gene_affected, transcript = data$transcript_affected, 
                           consequence = data$consequence_type, aa.change = data$aa_mutation)
  # save driver annotations only
  annotation <- annotation[which(annotation$gene %in% genes$ensembl.id), ]
  annotation$hugo.gene <- genes$gene_name[match(annotation$gene, genes$ensembl.id)]
  
  
  
  # assign contexts 
  vr <- VRanges(seqnames = Rle(paste('chr',data$chromosome, sep = '')), ranges = IRanges(start = data$chromosome_start, end = data$chromosome_end), 
                ref = as.character(data$mutated_from_allele), alt = as.character(data$mutated_to_allele), sampleNames = as.character(data$icgc_donor_id),
                submitted_sample_id = data$submitted_sample_id, strategy = data$sequencing_strategy)
  file = "/Users/[user]/data/ref/ucsc.hg19.fasta"
  fa <- open(FaFile(file, sprintf("%s.fai", file)))

  # check alignment genome is hg19
  vr.ref = DNAStringSet(ref(vr))
  file = "/Users/[user]/data/ref/ucsc.hg19.fasta"
  fa <- open(FaFile(file, sprintf("%s.fai", file)))
  gr <- granges(vr)
  genome.ref = getSeq(fa, gr)
  if (length(which(vr.ref != genome.ref)) > (length(vr) / 100)){
    error.names <- c(error.names, short.name)
    icgc.errors <- c(icgc.errors, 'non-hg19 alignment')
    next
  }
  
  vr <- mutationContext(vr = vr, ref = fa)
  vr$signature = as.factor(paste(vr$alteration, vr$context))
  vr$fullMutation <- paste(seqnames(vr), start(ranges(vr)), vr$signature)
  #vr$fullMutation <- paste(seqnames(vr), start(ranges(vr)), vr$signature, vr$aa_mutation)
  
  assign(paste(gsub('-', '_', short.name), '.ICGC.vr', sep = ''), vr)
  assign(paste(gsub('-', '_', short.name), '.ICGC.annotation', sep = ''), annotation)
  save(list = c(paste(gsub("-", "_", short.name), '.ICGC.vr', sep = ''), paste(gsub('-', '_', short.name), '.ICGC.annotation', sep = '')),
       file = paste('/Users/[user]/data/driver_spectrum/R_data/', short.name, '.ICGC.Rfile', sep = ''))
  rm(list = c(paste(gsub("-", "_", short.name), '.ICGC.vr', sep = ''), paste(gsub('-', '_', short.name), '.ICGC.annotation', sep = '')))
}

####################################################################################################################################
# add TCGA data 

tcga.errors <- c()
error.files <- c()
for (short.name in short.names){
  if (msl$TCGA.portal[match(short.name, msl$Study.ID)] == 'no'){
    next
  }
  print(short.name)
  
  tcga.abbr <- msl$TCGA.Study.Abbr[match(short.name, msl$Study.ID)]
  paths <- tcga.samples$Path[which(tcga.samples$Data.ID == tcga.abbr)]
  
  vr <- VRanges()
  annotation <- c()
  for (i in 1:length(paths)){
    file <- paste("/Users/[user]/apocrita/data_mount_point_2/daniel_import/signatures/TCGA/", paths[i], sep = '')
    
    if(class(try(read.table(file = file, sep = '\t',header = T, quote = '', stringsAsFactors = F), silent = T)) == 'try-error'){
      error.files <- c(error.files, file)
      tcga.errors <- c(tcga.errors, 'could not read')
      next
    }
    data <- read.table(file = file, sep = '\t',header = T, quote = '', stringsAsFactors = F)
    
    # SNP's only
    data <- data[which(data$Variant_Type == 'SNP'), ]
    
    # remove mitochondrial mutations
    data <- data[which(data$Chromosome %in% c(1:22, 'X', 'Y')), ]
    
    if(nrow(data) == 0){
      next
    }
    
    bases <- c('A', 'C', 'G', 'T')
    refs <- lapply(1:nrow(data), function(x) 
      c(data$Match_Norm_Seq_Allele1[x], data$Match_Norm_Seq_Allele2[x])[which(c(data$Match_Norm_Seq_Allele1[x], data$Match_Norm_Seq_Allele2[x]) %in% bases)])
    alts <- lapply(1:nrow(data), function(x) 
      c(data$Tumor_Seq_Allele1[x], data$Tumor_Seq_Allele2[x])[which(c(data$Tumor_Seq_Allele1[x], data$Tumor_Seq_Allele2[x]) %in% bases)])
    
    # all mutations
    mutation <- lapply(1:nrow(data), function(x) alts[[x]][which(!(alts[[x]] %in% refs[[x]]))])
    
    # remove cases with no mutation
    no.mutation <- which(sapply(mutation, length) == 0)
    if (length(no.mutation) > 0){
      data <- data[-no.mutation, ]
      mutation <- mutation[-no.mutation]
    }

    # check for multiple mutations in same row
    if (length(which(sapply(mutation, length) > 1)) > 0){
      error.files <- c(error.files, file)
      tcga.errors <- c(tcga.errors, 'Multiple Mutations in One Row')
    }
    
    # check for malformed bases in ref or alt
    if ((length(grep('(T|C|G|A)', refs)) != nrow(data)) | (length(grep('(T|C|G|A)', mutation)) != nrow(data))) {
      error.files <- c(error.files, file)
      tcga.errors <- c(tcga.errors, 'Malformed Ref Bases')
      next
    }
    
    # create annotations
    annotation.tmp <- data.frame(chromosome = data$Chromosome,  pos = data[, grep('(Start_Position|Start_position)' ,colnames(data))], 
                                 ref = sapply(refs, function(x) x[[1]][1]), alt = sapply(mutation, function(x) x[[1]][1]),
                                 dbSNP_RS = data$dbSNP_RS, consequence = data$Variant_Classification, hugo.gene = data$Hugo_Symbol)
    # save driver annotations only
    annotation.tmp <- annotation.tmp[which(annotation.tmp$hugo.gene %in% genes$gene_name), ]
    
    # cannot proceed if no data
    if(nrow(data) == 0){
      next
    }
    
    # find patient id's
    data$patient.id <- sapply(1:nrow(data), function(x) paste(strsplit(data$Tumor_Sample_Barcode[x], '-')[[1]][1:3], collapse = '-'))
    
    vr.tmp <- VRanges(seqnames = Rle(paste('chr',data$Chromosome, sep = '')), 
                      ranges = IRanges(start = data[, grep('(Start_Position|Start_position)' ,colnames(data))], end = data[, grep('(End_Position|End_position)' ,colnames(data))]), 
                      ref = as.character(data$Match_Norm_Seq_Allele1), alt = as.character(data$Tumor_Seq_Allele2), 
                      sampleNames = as.character(data$patient.id), submitted_sample_id = data$Tumor_Sample_Barcode, strategy = data$Sequence_Source)
  
    
    # check alignment genome is hg19
    vr.tmp.ref = DNAStringSet(ref(vr.tmp))
    file = "/Users/[user]/data/ref/ucsc.hg19.fasta"
    fa <- open(FaFile(file, sprintf("%s.fai", file)))
    gr <- granges(vr.tmp)
    genome.ref = getSeq(fa, gr)
    if (length(which(vr.tmp.ref != genome.ref)) > (length(vr.tmp) / 100)){
      error.files <- c(error.files, paths[i])
      tcga.errors <- c(tcga.errors, 'non-hg19 alignment')
      next
    }
    
    vr <- c(vr, vr.tmp)
    annotation <- rbind(annotation, annotation.tmp)
  }
  if(length(vr) == 0){
    next
  }
  # remove duplicates from the vr and annotation file
  
  # annotate contexts
  file = "/Users/[user]/data/ref/ucsc.hg19.fasta"
  fa <- open(FaFile(file, sprintf("%s.fai", file)))
  vr <- mutationContext(vr = vr, ref = fa)
  vr$signature = as.factor(paste(vr$alteration, vr$context))
  vr$fullMutation <- paste(seqnames(vr), start(ranges(vr)), vr$signature)
  
  assign(paste(gsub('-', '_', short.name), '.TCGA.vr', sep = ''), vr)
  assign(paste(gsub('-', '_', short.name), '.TCGA.annotation', sep = ''), annotation)
  save(list = c(paste(gsub("-", "_", short.name), '.TCGA.vr', sep = ''), paste(gsub('-', '_', short.name), '.TCGA.annotation', sep = '')),
       file = paste('/Users/[user]/data/driver_spectrum/R_data/', short.name, '.TCGA.Rfile', sep = ''))
  rm(list = c(paste(gsub("-", "_", short.name), '.TCGA.vr', sep = ''), paste(gsub('-', '_', short.name), '.TCGA.annotation', sep = '')))
}

####################################################################################################################################
# merge the data sources and add amino acids to annotation
for (short.name in short.names){
  print(paste('Mergeing ', short.name, '...', sep = ''))
  icgc.file <- paste('/Users/[user]/data/driver_spectrum/R_data/', short.name, '.ICGC.Rfile', sep = '')
  tcga.file <- paste('/Users/[user]/data/driver_spectrum/R_data/', short.name, '.TCGA.Rfile', sep = '')
  merged.file <- paste('/Users/[user]/data/driver_spectrum/R_data/', short.name, '.merged.Rfile', sep = '')
  if (file.exists(icgc.file) & !file.exists(tcga.file)){
    load(file = icgc.file)
    icgc.annotation <- get(paste(gsub('-', '_', short.name), '.ICGC.annotation', sep = ''))
    merged.annotation <- icgc.annotation[ , c('chromosome', 'pos', 'ref', 'alt', 'consequence')]
    icgc.vr <- get(paste(gsub('-', '_', short.name), '.ICGC.vr', sep = ''))
    merged.vr <- icgc.vr
    rm(list = c(paste(gsub("-", "_", short.name), '.ICGC.vr', sep = ''), paste(gsub('-', '_', short.name), '.ICGC.annotation', sep = '')))
  } else if (file.exists(tcga.file) & !file.exists(icgc.file)){
    load(file = tcga.file)
    tcga.annotation <- get(paste(gsub('-', '_', short.name), '.TCGA.annotation', sep = ''))
    merged.annotation <- tcga.annotation[ , c('chromosome', 'pos', 'ref', 'alt', 'consequence')]
    tcga.vr <- get(paste(gsub('-', '_', short.name), '.TCGA.vr', sep = ''))
    merged.vr <- tcga.vr
    rm(list = c(paste(gsub("-", "_", short.name), '.TCGA.vr', sep = ''), paste(gsub('-', '_', short.name), '.TCGA.annotation', sep = '')))
  } else if (file.exists(tcga.file) & file.exists(icgc.file)){
    load(file = icgc.file)
    load(file = tcga.file)
    icgc.annotation <- get(paste(gsub('-', '_', short.name), '.ICGC.annotation', sep = ''))
    tcga.annotation <- get(paste(gsub('-', '_', short.name), '.TCGA.annotation', sep = ''))
    
    # create merged.annotation
    merged.annotation <- rbind(icgc.annotation[ , c('chromosome', 'pos', 'ref', 'alt', 'consequence')],
                               tcga.annotation[ , c('chromosome', 'pos', 'ref', 'alt', 'consequence')])
    
    # merge mutation tables
    icgc.vr <- get(paste(gsub('-', '_', short.name), '.ICGC.vr', sep = ''))
    tcga.vr <- get(paste(gsub('-', '_', short.name), '.TCGA.vr', sep = ''))
    
    tcga.samples <- unique(as.character(sampleNames(tcga.vr)))
    new.samples <- tcga.samples[which(!sapply(1:length(tcga.samples), function(x) length(grep(tcga.samples[x], icgc.vr$submitted_sample_id))) > 0)]
    
    levels(sampleNames(icgc.vr)) <- c(levels(sampleNames(icgc.vr)), new.samples)
    
    merged.vr <- c(icgc.vr, tcga.vr[which(sampleNames(tcga.vr) %in% new.samples), ])
    rm(list = c(paste(gsub("-", "_", short.name), '.ICGC.vr', sep = ''), paste(gsub('-', '_', short.name), '.ICGC.annotation', sep = '')))
    rm(list = c(paste(gsub("-", "_", short.name), '.TCGA.vr', sep = ''), paste(gsub('-', '_', short.name), '.TCGA.annotation', sep = '')))
  } else{
    next
  }
  
  # keep mutations from most common strategy for each sample only
  tmp <- table(sampleNames(merged.vr), merged.vr$strategy)
  informative.strategy <- data.frame(sample = levels(sampleNames(merged.vr)))
  informative.strategy$strategy <- apply(tmp, 1, function(x) colnames(tmp)[which.max(x)])
  merged.vr$informative.strategy <- informative.strategy$strategy[match(sampleNames(merged.vr), informative.strategy$sample)]
  merged.vr <- merged.vr[which(merged.vr$strategy == merged.vr$informative.strategy),]

  # remove duplicates from merged vr
  mutation <- paste(sampleNames(merged.vr), merged.vr$fullMutation)
  merged.vr <- merged.vr[match(unique(mutation), mutation), ]
  
  # annotate the merged annotation file
  # remove duplicates - not done with other mutations
  merged.annotation$mutation <- paste(merged.annotation$chromosome, merged.annotation$pos, merged.annotation$ref, merged.annotation$alt)
  merged.annotation <- merged.annotation[match(unique(merged.annotation$mutation), merged.annotation$mutation), ]
  # annotate A.A. changes 
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  file = "/Users/[user]/data/ref/ucsc.hg19.fasta"
  fa <- open(FaFile(file, sprintf("%s.fai", file)))
  merged.annotation.vr <- VRanges(seqnames = Rle(paste('chr', merged.annotation$chromosome, sep = '')), 
                                ranges = IRanges(start = merged.annotation$pos, end = merged.annotation$pos), 
                                ref = merged.annotation$ref, alt = merged.annotation$alt, mutation = merged.annotation$mutation)
  vr <- mutationContext(vr = merged.annotation.vr, ref = fa)
  vr$signature = as.factor(paste(vr$alteration, vr$context))
  vr$fullMutation <- paste(seqnames(vr), start(ranges(vr)), vr$signature)
  # copy over full mutation
  merged.annotation$fullMutation <- vr$fullMutation[match(merged.annotation$mutation, vr$mutation)]
  # remove temporary id column
  merged.annotation$mutation <- NULL
  assign(paste(gsub('-', '_', short.name), '.merged.annotation', sep = ''), merged.annotation)
  
  assign(paste(gsub('-', '_', short.name), '.merged.vr', sep = ''), merged.vr)
  
  save(list = c(paste(gsub("-", "_", short.name), '.merged.vr', sep = ''), paste(gsub('-', '_', short.name), '.merged.annotation', sep = '')),
       file = merged.file)
  rm(list = c(paste(gsub("-", "_", short.name), '.merged.vr', sep = ''), paste(gsub("-", "_", short.name), '.merged.annotation', sep = '')))
}

####################################################################################################################################

# find the oncogene hits and frequent oncogene hits in each set [VERY SLOW]
# create inputs for annovar
hits.all <- data.frame()
for (short.name in short.names){
  file = paste('/Users/[user]/data/driver_spectrum/R_data/', short.name, '.merged.Rfile', sep = '')
  if (!file.exists(file)){
    next
  }
  load(file)
  print(paste('Calculating', short.name))
  hits.all.tmp <- get(paste(gsub('-', '_', short.name), 'merged.annotation', sep = '.'))
  hits.all <- rbind(hits.all, hits.all.tmp)
  rm(list = c(paste(gsub("-", "_", short.name), '.merged.vr', sep = ''), paste(gsub("-", "_", short.name), '.merged.annotation', sep = '')))
}
# remove duplicate rows from hits.all
hits.all <- hits.all[match(unique(hits.all$fullMutation), hits.all$fullMutation),]
# write annovar input
annovar.input <- data.frame(CHROM = hits.all$chromosome, POS = hits.all$pos, POS = hits.all$pos, ref = hits.all$ref, alt = hits.all$alt, fullMutation = hits.all$fullMutation)
write.table(annovar.input, file = '/Users/[user]/data/driver_spectrum/driver_mutations/annovar.csv', quote = F, row.names = F, sep = '\t')

# read in the annotated file and update hits.all
df <- read.table(file = '/Users/[user]/data/driver_spectrum/driver_mutations/annovar.csv.exonic_variant_function', sep = '\t')
hits.all <- data.frame(chromosome = df[,4], pos = df[,5], ref = df[,7], alt = df[,8], consequence = df[,2], fullMutation = df[,9], annovar = df[,3])
hits.all$annovar <- as.character(hits.all$annovar)
hits.all$transcripts <- sapply(1:nrow(hits.all), function(x) ParseAnnovarTranscript(string = hits.all$annovar[x]))
hits.all$mutations <- sapply(1:nrow(hits.all), function(x) ParseAnnovarMutation(string = hits.all$annovar[x]))
hits.all$hugo.gene <- sapply(hits.all$annovar, function(x) ParseAnnovarGene(x))
hits.all$annovar <- NULL
# driver genes only
hits.all <- hits.all[which(hits.all$hugo.gene %in% genes$gene_name), ]
# non-silent only 
hits.all <- hits.all[which(hits.all$consequence %in% c('nonsynonymous SNV', 'stopgain', 'stoploss')), ]
hits.all$transcript <- genes$transcript[match(hits.all$hugo.gene, genes$gene_name)]
hits.all$residue <- sapply(1:nrow(hits.all), function(x) strsplit(hits.all$mutations[x], ":")[[1]][match(hits.all$transcript[x], strsplit(hits.all$transcripts[x], ':')[[1]])])
# remove mutations with no effect on canonical transcript
hits.all <- hits.all[which(!is.na(hits.all$residue)), ]
hits.all$name <- paste(hits.all$hugo.gene, hits.all$residue)
hits.all$fullMutation <- as.character(hits.all$fullMutation)
hits.all$context <- sapply(hits.all$fullMutation, function(x) paste(strsplit(x, ' ')[[1]][3:4], collapse = ' '))
hits.all$gene.context <- paste(hits.all$hugo.gene, hits.all$context)
hits.all$aa.context <- paste(hits.all$name, ' (', hits.all$context, ')', sep = '')
hits.all$classification <- genes$classification[match(hits.all$hugo.gene, genes$gene_name)]
hits.all$context.name <- sapply(1:nrow(hits.all), function(x) paste(sapply(strsplit(hits.all$context[x], ';')[[1]], function(y) paste(substr(y,4,4), substr(y,1,1), substr(y,6,6), '>', substr(y,4,4), substr(y,2,2), substr(y,6,6), sep = '')), collapse = ':'))
hits.all$short.name <- paste(hits.all$name, ' (', hits.all$context.name, ')', sep = '')

save(list = c('hits.all'),
     file = paste('/Users/[user]/data/driver_spectrum/R_data/', 'Hits.Rfile', sep = ''))

####################################################################################################################################

# calculate signature and assess hits for each patient
signature.table.all <- c()
patient.metadata <- data.frame()
presence.table <- c()
for (short.name in short.names){
  print(paste('Calculating', short.name))
  file = paste('/Users/[user]/data/driver_spectrum/R_data/', short.name, '.merged.Rfile', sep = '')
  if (!file.exists(file)){
    next
  }
  load(file = paste('/Users/[user]/data/driver_spectrum/R_data/', short.name, '.merged.Rfile', sep = ''))
  vr <- get(paste(gsub('-', '_', short.name), '.merged.vr', sep = ''))
  sampleNames(vr) <- factor(sampleNames(vr), levels = unique(sampleNames(vr)))
  # create signature.table.tmp
  
  vr$signature <- factor(vr$signature, levels = processes)
  signature.table.tmp <- table(vr$signature, sampleNames(vr))
  patient.metadata.tmp <- data.frame(patient.ID = colnames(signature.table.tmp), 
                                     short.name = short.name, 
                                     disease = msl$Disease[match(short.name, msl$Study.ID)],
                                     disease.abbr = msl$Disease.Abbr[match(short.name, msl$Study.ID)],
                                     strategy = vr$strategy[match(colnames(signature.table.tmp), sampleNames(vr))],
                                     tcga.id = NA)
  tcga.ids <- grep('TCGA', vr$submitted_sample_id[match(rownames(patient.metadata.tmp), sampleNames(vr))])
  if(length(tcga.ids) > 0){
    patient.metadata.tmp$tcga.id[tcga.ids] <- sapply(vr$submitted_sample_id[match(rownames(patient.metadata.tmp), sampleNames(vr))[tcga.ids]], function(x) paste(strsplit(x, '-')[[1]][1:3], collapse = '-'))
  }
  # add to singature.table and signature.table.metadata
  signature.table.all <- cbind(signature.table.all, signature.table.tmp)
  patient.metadata <- rbind(patient.metadata, patient.metadata.tmp)

  # create presence.table.tmp
  sub <- vr[which(vr$fullMutation %in% hits.all$fullMutation), ]
  sub$fullMutation <- factor(sub$fullMutation, levels = unique(hits.all$fullMutation))
  presence.table.tmp <- table(sub$fullMutation, sampleNames(sub))
  presence.table <- cbind(presence.table, presence.table.tmp)
  rm(list = c(paste(gsub("-", "_", short.name), '.merged.vr', sep = ''), paste(gsub("-", "_", short.name), '.merged.annotation', sep = '')))
}
presence.table <- as.data.frame(presence.table)
patient.metadata$disease.group <- cancer.groupings$disease.group[match(patient.metadata$disease, cancer.groupings$Disease)]
save(list = c('signature.table.all', 'presence.table', 'patient.metadata'),
     file = paste('/Users/[user]/data/driver_spectrum/R_data/', 'ByPatientSummaries.Rfile', sep = ''))

####################################################################################################################################

# create gene.presence.table and gene.context.presence.table
tmp <- as.data.frame(which(presence.table != 0, arr.ind = T))
tmp$patient <- factor(colnames(presence.table)[tmp$col], levels = colnames(presence.table))
tmp$gene <- hits.all$hugo.gene[match(rownames(tmp), hits.all$fullMutation)]
tmp$name <- hits.all$name[match(rownames(tmp), hits.all$fullMutation)]
#tmp$gene <- sapply(rownames(tmp), function(x) strsplit(x, ' ')[[1]][1])
tmp$context <- sapply(rownames(tmp), function(x) paste(strsplit(x, ' ')[[1]][3:4], collapse = ' '))
tmp$gene.context <- paste(tmp$gene, tmp$context)
tmp$aa.context <- hits.all$aa.context[match(rownames(tmp), hits.all$fullMutation)]
context.presence.table <- table(tmp$context, tmp$patient)
gene.context.presence.table <- table(tmp$gene.context, tmp$patient)
gene.presence.table <- table(tmp$gene, tmp$patient)
aa.presence.table <- table(tmp$name, tmp$patient)
aa.context.presence.table <- table(tmp$aa.context, tmp$patient)

signature.table <- signature.table.all - context.presence.table
patient.metadata$signature.mutations <- colSums(signature.table)
mutations.and.patients <- tmp

save(list = c('patient.metadata', 'signature.table', 'mutations.and.patients', 'gene.context.presence.table', 'gene.presence.table', 'aa.presence.table', 'aa.context.presence.table'),
     file = paste('/Users/[user]/data/driver_spectrum/R_data/', 'Presence.Rfile', sep = ''))
