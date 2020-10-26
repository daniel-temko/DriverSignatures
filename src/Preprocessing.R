require(reshape)
require(reshape2)
require(nnls)
require(RColorBrewer)
require(Hmisc)
require(fpc)

AggregateOccurencesByCancer <- function(table, metadata, group.by){
  data <- data.frame(row.names = row.names(table))
  ids <- unique(metadata[[group.by]])
  for (id in ids){
    sub <- table[ , match(metadata$patient.ID[which(metadata[[group.by]] == id)], colnames(table))]
    if (length(which(metadata[[group.by]] == id)) > 1) {
      sums <- rowSums(sub > 0)
    } else if (length(which(metadata[[group.by]] == id)) == 1){
      sums <- sub > 0
    } else {
      sums <- rep(0, nrow(table))
    }
    data <- cbind(data, sums)
  }
  colnames(data) <- ids
  return(data)
}

AggregateByCancer <- function(table, metadata, group.by){
  data <- data.frame(row.names = row.names(table))
  ids <- levels(metadata[[group.by]])
  for (id in ids){
    sub <- table[ , match(metadata$patient.ID[which(metadata[[group.by]] == id)], colnames(table))]
    if (length(which(metadata[[group.by]] == id)) > 1) {
      sums <- rowSums(sub)
    } else if (length(which(metadata[[group.by]] == id)) == 1){
      sums <- sub
    } else {
      sums <- rep(0, nrow(table))
    }
    data <- cbind(data, sums)
  }
  colnames(data) <- ids
  return(data)
}

RestrictedRegression <- function(mutations, signature.map, signature.group, strategy, prc.proven.exome, prc.proven){
  if (signature.group %in% signature.map$signature.group){
    sub <- signature.map[which(signature.map$signature.group == signature.group), ]
    indices <- which(colnames(prc.proven) %in% sub$variable)
    if(strategy == 'WXS'){
      reg <- nnls(as.matrix(prc.proven.exome[, indices]), mutations)
    } else if(strategy == 'WGS'){
      reg <- nnls(as.matrix(prc.proven[, indices]), mutations)
    } else {
      warning('Strategy is not \'WXS\' or \'WGS\'')
    }
    res <- rep(0, ncol(prc.proven))
    res[indices] <- reg$x
    return(res)
  } else {
    return(rep(0, ncol(prc.proven)))
  }
}


ListMutations <- function(genes, threshold){
  hits.sub <- hits.all[which((hits.all$hugo.gene %in% genes) & (hits.all$aa.context %in% rownames(mutation.aa.context.tab))), ]
  hits.sub <- hits.sub[unique(match(hits.sub$aa.context, hits.sub$aa.context)),]
  
  full.list <- hits.sub$aa.context
  counts <- rowSums(mutation.aa.context.tab[match(full.list, rownames(mutation.aa.context.tab)),])
  list <- full.list[which((counts / sum(counts)) >= threshold)]
  return(list)
}


DifferentialSelection <- function(ids, mutations){
  # create a 96-channel signature based on even tnf's
  sig.scores <- signature.scores.norm[, match(ids, colnames(signature.scores.norm))]
  sums <- rowSums(sig.scores)
  sums <- sums / sum(sums)
  tmp <- t(as.matrix(prc.proven.even) %*% as.matrix(sums))
  sds <- apply(t(as.matrix(prc.proven.even) %*% as.matrix(sig.scores)), 2, sd)
  
  # find all non-silent mutations in the gene
  hits.sub <- hits.all[match(mutations, hits.all$aa.context), ]
  
  # create a data frame of pm and nm for each mutation
  res <- data.frame(mut = mutations, context = hits.sub$context, p = tmp[match(hits.sub$context, processes)], sd = sds[match(hits.sub$context, processes)])
  mut.sub <- mutation.aa.context.tab[match(res$mut, rownames(mutation.aa.context.tab)), match(ids, colnames(mutation.aa.context.tab))]
  res$n <- rowSums(mut.sub)[match(res$mut, rownames(mut.sub))]
  
  res$enrichment <- res$n / res$p
  res$normalised.enrichment <- res$enrichment / sum(res$enrichment)
  return(res)
}

TransformContext <- function(context){
  paste(substr(context,4,4), substr(context,1,1), substr(context,6,6), '>', substr(context,4,4), substr(context,2,2), substr(context,6,6), sep = '')
}

# signature metadata - used for assigning signature scores and labelling signatures
signature.annotations <- read.table(file = "signature.annotations.csv", 
                                    sep = ',', col.names = c('process', 'association', 'type', 'clock.like'))
signature.annotations$name <- paste(signature.annotations$process, signature.annotations$association)
signature.annotations$unq.ann <- as.character(signature.annotations$association)
signature.annotations$unq.ann[which(signature.annotations$unq.ann == 'Unknown aetiology')] <- signature.annotations$name[which(signature.annotations$unq.ann == 'Unknown aetiology')]
signatures.tissue.map <- read.table(file = "diseases.and.signatures.csv", 
                                    sep = ',', header = T)
signature.tissue.df <- melt(signatures.tissue.map, id.var = 'signature.group')
signature.tissue.df <- signature.tissue.df[-which(is.na(signature.tissue.df$value)), ]
signature.tissue.df$unq.ann <- signature.annotations$unq.ann[match(signature.tissue.df$variable, signature.annotations$process)]
signature.tissue.df$value <- NULL

min.signature.mutations <- 20
freq.threshold <- 0.01
n.base <- 10

load(file = paste(wrk_dir, 'R_data/signatures.summary.data.Rfile', sep = ''))

####################################################################################################################################

## filter data
diseases.with.signatures <- cancer.groupings$Disease[which(!is.na(cancer.groupings$Signature.Group))]
init.ids <- patient.metadata$patient.ID[which((patient.metadata$strategy %in% c('WGS', 'WXS')) &
                                               (msl$Disease[match(patient.metadata$short.name, msl$Study.ID)] %in% diseases.with.signatures))]


temp.ids <- intersect(init.ids, patient.metadata$patient.ID[which((patient.metadata$signature.mutations >= min.signature.mutations))])
metadata <- patient.metadata[which(patient.metadata$patient.ID %in% temp.ids), ]
signature.context.tab <- as.data.frame.matrix(signature.table[ ,  which(colnames(signature.table) %in% temp.ids)])
signature.context.norm <- apply(signature.context.tab, 2, function(x) x / sum(x))


metadata <- patient.metadata[which(patient.metadata$patient.ID %in% temp.ids), ]
totals <- table(metadata$disease.group)

min.cases <- 100
eligible.cancers <- names(totals)[which(totals > min.cases )]
ids.full <- intersect(init.ids, patient.metadata$patient.ID[which(patient.metadata$disease.group %in% eligible.cancers)])
ids <- intersect(ids.full, temp.ids)

# restrict to patients with sufficient mutations for signature assingment
metadata.full <- patient.metadata[which(patient.metadata$patient.ID %in% ids.full), ]
mutation.tab.full <- presence.table[ , which(colnames(presence.table) %in% ids.full)]
mutation.context.tab.full <- gene.context.presence.table[ , which(colnames(gene.context.presence.table) %in% ids.full)]
mutation.gene.tab.full <- gene.presence.table[ , which(colnames(gene.presence.table) %in% ids.full)]
mutation.aa.change.tab.full <- aa.presence.table[ , which(colnames(aa.presence.table) %in% ids.full)]
mutation.aa.context.tab.full <- aa.context.presence.table[ , which(colnames(aa.presence.table) %in% ids.full)]
signature.context.tab.full <- signature.table[ ,  which(colnames(signature.table) %in% ids.full)]


metadata <- patient.metadata[which(patient.metadata$patient.ID %in% ids), ]
mutation.tab <- presence.table[ , which(colnames(presence.table) %in% ids)]
mutation.context.tab <- gene.context.presence.table[ , which(colnames(gene.context.presence.table) %in% ids)]
mutation.gene.tab <- gene.presence.table[ , which(colnames(gene.presence.table) %in% ids)]
mutation.aa.change.tab <- aa.presence.table[ , which(colnames(aa.presence.table) %in% ids)]
mutation.aa.context.tab <- aa.context.presence.table[ , which(colnames(aa.presence.table) %in% ids)]
signature.context.tab <- as.data.frame.matrix(signature.table[ ,  which(colnames(signature.table) %in% ids)])
signature.context.norm <- apply(signature.context.tab, 2, function(x) x / sum(x))
strat.table <- data.frame(strategy = c('WXS', 'WGS'), name = c('exome', 'genome'))
metadata$disease.group <- factor(metadata$disease.group)

# use nnls to find signature scores
signature.scores.tab <- as.data.frame(lapply(1:ncol(signature.context.tab), function(x) 
  RestrictedRegression(mutations = signature.context.tab[,x],
                       signature.map = signature.tissue.df,
                       signature.group = cancer.groupings$Signature.Group[match(as.character(metadata$disease[match(colnames(signature.context.tab)[x], metadata$patient.ID)]), cancer.groupings$Disease)],
                       strategy = metadata$strategy[match(colnames(signature.context.tab)[x], metadata$patient.ID)],
                       prc.proven.exome = prc.proven.exome,
                       prc.proven = prc.proven)))
signature.rates <- as.data.frame(lapply(1:nrow(metadata), function(x) if(metadata$strategy[x] == 'WGS') {signature.scores.tab[,x] / genome.signature.scaling.factors} else {signature.scores.tab[,x] / exome.signature.scaling.factors}))
colnames(signature.rates) <- metadata$patient.ID
rownames(signature.rates) <- paste('Signature.', 1:30, sep = '')
signature.scores.norm <- apply(signature.rates, 2, function(x) x / sum(x) )

contexts.imputed <- as.matrix(prc.proven.even) %*% as.matrix(signature.scores.norm)

selection.genes <- c('TP53', 'KRAS', 'NRAS', 'PIK3CA', 'CTNNB1', 'SMAD4', 'BRAF', 'IDH1', 'IDH2', 'gs.1', 'gs.2', 'gs.3') 
gene.sets <- data.frame(id = paste('gs.', 1:3, sep = ''), 
                        gene.1 = c('KRAS', 'CTNNB1', 'IDH1'), 
                        gene.2 = c('BRAF', 'APC', 'IDH2'), 
                        gene.3 = c('NRAS', NA, NA),
                        name = c('KRAS, BRAF, and NRAS', 'APC and CTNNB1', 'IDH1 and IDH2'))

data <- AggregateOccurencesByCancer(mutation.aa.context.tab, metadata, 'disease.group')
set.diseases <- c()
for(gene in selection.genes){
  print(gene)
  if(gene %in% gene.sets$id){
    set.genes <- gene.sets[match(gene, gene.sets$id),-1][which(!is.na(gene.sets[match(gene, gene.sets$id),-1]))]
    mutations <- unname(unlist(sapply(set.genes, function(x) ListMutations(genes = x, threshold = freq.threshold))))
  } else {
    mutations <- ListMutations(gene, freq.threshold)
  }
  tmp <- data.frame(disease = levels(metadata$disease.group), set = gene)
  tmp$counts <- apply(data[match(mutations, rownames(data)), match(tmp$disease, colnames(data))], 2, sum)
  set.diseases <- rbind(set.diseases, tmp)
}
set.diseases <- set.diseases[which(set.diseases$counts > n.base), ]

#save
save(list = c('metadata.full', 'mutation.aa.change.tab', 'mutation.aa.context.tab', 'signature.scores.tab', 'signature.rates', 'signature.scores.norm',
              'signature.context.tab', 'signature.context.norm', 'contexts.imputed', 'signature.tissue.df',
              'hits.all', 'metadata', 'selection.genes', 'gene.sets', 'set.diseases', 'processes',
              'prc.proven', 'prc.proven.even', 'prc.proven.exome', 'freq.threshold'),
     file = paste(wrk_dir, 'R_data/signatures.summary.data.Rfile', sep = ''))
