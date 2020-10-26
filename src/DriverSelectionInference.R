require(ggrepel)
require(poibin)

TestForDifference <- function(mut1, mut2, mut.sub){
  mut1.process <- hits.all$context[match(mut1, hits.all$aa.context)]
  mut2.process <- hits.all$context[match(mut2, hits.all$aa.context)]
  
  mut1.ids <- colnames(mut.sub)[which(mut.sub[match(mut1, rownames(mut.sub)),] == 1)]
  mut2.ids <- colnames(mut.sub)[which(mut.sub[match(mut2, rownames(mut.sub)),] == 1)]
  mut.ids <- c(mut1.ids, mut2.ids)
  if(length(mut.ids) == 0){return(rep(NA, 3))}
  
  mut1.prob <- contexts.imputed[match(mut1.process, processes), match(mut.ids, colnames(contexts.imputed))]
  mut2.prob <- contexts.imputed[match(mut2.process, processes), match(mut.ids, colnames(contexts.imputed))]
  
  mut1.n <- length(mut1.ids)
  
  threshold.prob <- dpoibin(kk = mut1.n, pp = mut1.prob / (mut1.prob + mut2.prob))
  probs <- dpoibin(kk = 0:length(mut1.prob), pp = mut1.prob / (mut1.prob + mut2.prob))
  p.value <- sum(probs[which(probs <= threshold.prob)])
  p.lower <- sum(dpoibin(kk = 0:mut1.n, pp = mut1.prob / (mut1.prob + mut2.prob)))
  p.higher <- sum(dpoibin(kk = mut1.n:length(mut1.prob), pp = mut1.prob / (mut1.prob + mut2.prob)))
  
  return(c(p.value, p.lower, p.higher))
}

RelativeRisk <- function(mut1, mut2, mut.sub, bootstrap.iterations){
  mut1.process <- hits.all$context[match(mut1, hits.all$aa.context)]
  mut2.process <- hits.all$context[match(mut2, hits.all$aa.context)]
  
  mut1.ids <- colnames(mut.sub)[which(mut.sub[match(mut1, rownames(mut.sub)),] == 1)]
  mut2.ids <- colnames(mut.sub)[which(mut.sub[match(mut2, rownames(mut.sub)),] == 1)]
  mut.ids <- c(mut1.ids, mut2.ids)
  if(length(mut.ids) == 0){return(rep(NA, 3))}
  
  mut1.prob <- contexts.imputed[match(mut1.process, processes), match(mut.ids, colnames(contexts.imputed))]
  mut2.prob <- contexts.imputed[match(mut2.process, processes), match(mut.ids, colnames(contexts.imputed))]
  
  mut1.n <- length(mut1.ids)
  
  # relative risk estimate
  rr.est <- EstimateRR(mut1.n, mut1.prob, mut2.prob)
  if(is.na(rr.est) | is.infinite(rr.est)){return(c(rr.est, rep(NA, 2)))}
  
  # confidence interval
  pp.est <- ((mut1.prob * rr.est) / ((mut1.prob * rr.est) + mut2.prob))
  mut1.n.sim <- rpoibin(bootstrap.iterations, pp = pp.est)
  rr.simulated <- sapply(mut1.n.sim, function(x) EstimateRR(x, mut1.prob, mut2.prob))
  lb <- as.numeric(quantile(rr.simulated, 0.025, na.rm  = T))
  ub <- as.numeric(quantile(rr.simulated, 0.975, na.rm = T))
  return(c(rr.est, lb, ub))
}

EstimateRR <- function(mut1.n, mut1.prob, mut2.prob){
  if(mut1.n == 0) {return(0)}
  if(mut1.n == length(mut1.prob)) {return(Inf)}
  ub <- 1
  while((DLLDR(mut1.n, mut1.prob, mut2.prob, ub) > 0) | (DLLDR(mut1.n, mut1.prob, mut2.prob, 1e-20) < 0)){
    ub <- ub * 10
    if(ub > 10^10){return(NA)}
  }
  uni <- uniroot(function(x) DLLDR(mut1.n, mut1.prob, mut2.prob, x), c(1e-20, ub))$root
  return(uni)
}

DLLDR <- function(mut1.n, mut1.prob, mut2.prob, R){
  dlldr <- mut1.n/R - sum(mut1.prob / ((mut1.prob * R) + mut2.prob))
  return(dlldr)
}

AnalyseSelection <- function(mut.sub){
  df <- data.frame(mutation = rownames(mut.sub))
  df$process <- hits.all$context[match(df$mutation, hits.all$aa.context)]
  df$mean.p <- apply(contexts.imputed[match(df$process, processes), match(colnames(mut.sub), colnames(contexts.imputed))], 1, mean)
  df$n <- rowSums(mut.sub)
  return(df)
}

####################################################################################
# plotting parameters
width = 8.18
height = 4.9
all.cols <- c(rainbow(12)[c(11,7,3,2,1,5,9)])
pch.ids <- c(16, 15, 17, 18:25)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector[c(7, 13)] <- col_vector[26 + c(3, 4)]
color.map <- data.frame(prc = colnames(prc.proven), 
                        col = col_vector[match(signature.annotations$unq.ann, unique(signature.annotations$unq.ann))])
size.1 <- 2
size.2 <- 1.5

freq.threshold <- 0.01

####################################################################################

# pairwise analysis
pairwise.tests <- c()
for(gene in selection.genes){
  print(gene)
  if(gene %in% gene.sets$id){
    set.genes <- gene.sets[match(gene, gene.sets$id),-1][which(!is.na(gene.sets[match(gene, gene.sets$id),-1]))]
    mutations <- unname(unlist(sapply(set.genes, function(x) ListMutations(genes = x, threshold = freq.threshold))))
  } else {
    mutations <- ListMutations(gene, freq.threshold)
  }
  diseases <- unique(set.diseases$disease[which(set.diseases$set == gene)])
  for(disease in diseases){
    print(disease)
    res <- cbind(expand.grid(mutation1 = mutations, mutation2 = mutations), 
                 expand.grid(mutation1.id = 1:length(mutations), mutation2.id = 1:length(mutations), 
                             set = gene, 
                             disease = disease))
    ids <- metadata$patient.ID[which(metadata$disease.group == disease)]
    mut.sub <- mutation.aa.context.tab[match(mutations, rownames(mutation.aa.context.tab)),
                                       match(ids, colnames(mutation.aa.context.tab))]
    # exclude cases with more than one mutation from the set
    mut.sub <- mut.sub[, which(colSums(mut.sub) == 1)]
    res$gene1 <- hits.all$hugo.gene[match(res$mutation1, hits.all$aa.context)]
    res$gene2 <- hits.all$hugo.gene[match(res$mutation2, hits.all$aa.context)]
    res <- res[which(res$mutation1.id > res$mutation2.id),]
    
    tmp <- unname(t(as.data.frame(lapply(1:nrow(res), function(x)
      TestForDifference(mut1 = res$mutation1[x], 
                        mut2 = res$mutation2[x], 
                        mut.sub = mut.sub)))))
    res <- cbind(res, tmp)
    colnames(res)[(ncol(res) - 2):ncol(res)] <- c('p.value', 'p.lower', 'p.higher')
    res <- res[which(!is.na(res$p.lower)),]
    res$direction <- sapply(1:nrow(res), function(x) if(res$p.lower[x] < res$p.higher[x]) {'below'} else {'above'})
    
    res$n.mutation1 <- sapply(1:nrow(res), function(x) length(which((mut.sub[match(res$mutation1[x], rownames(mut.sub)),] == 1) & (mut.sub[match(res$mutation2[x], rownames(mut.sub)),] == 0))))
    res$n.mutation2 <- sapply(1:nrow(res), function(x) length(which((mut.sub[match(res$mutation2[x], rownames(mut.sub)),] == 1) & (mut.sub[match(res$mutation1[x], rownames(mut.sub)),] == 0))))
                              
    pairwise.tests <- rbind(pairwise.tests, res[,c('mutation1', 'n.mutation1', 'mutation2', 'n.mutation2', 'set', 'disease', 'p.value', 'direction')])
  }
}

fdr <- 0.05
selection.sum <- pairwise.tests

selection.sum <- selection.sum[order(selection.sum$p.value),]
selection.sum$stat <- sapply(1:nrow(selection.sum), function(x) (selection.sum$p.value[x] * length(which(!is.na(selection.sum$p.value)))) / length(which(selection.sum$p.value <= selection.sum$p.value[x])))
selection.sum$q.value <- sapply(1:nrow(selection.sum), function(x) if(is.na(selection.sum$stat[x])) {NA} else {min(selection.sum$stat[x:nrow(selection.sum)], na.rm = T)})
selection.sum$significant <- sapply(selection.sum$q.value, function(x) if(!is.na(x) & (x < fdr)) {'s'} else {'ns'})

selection.sum$gene1 <- hits.all$hugo.gene[match(selection.sum$mutation1, hits.all$aa.context)]
selection.sum$gene2 <- hits.all$hugo.gene[match(selection.sum$mutation2, hits.all$aa.context)]
save(list = c('selection.sum'),
     file = paste(wrk_dir, 'R_processed_data/', 'SelectionSum.Rfile', sep = ''))

selection.write <- selection.sum
selection.write$mutation1 <- sapply(1:nrow(selection.sum), function(x) if(selection.sum$direction[x] == 'above') {selection.sum$mutation1[x]} else {selection.sum$mutation2[x]})
selection.write$gene1 <- sapply(1:nrow(selection.sum), function(x) if(selection.sum$direction[x] == 'above') {selection.sum$gene1[x]} else {selection.sum$gene2[x]})
selection.write$mutation2 <- sapply(1:nrow(selection.sum), function(x) if(selection.sum$direction[x] == 'above') {selection.sum$mutation2[x]} else {selection.sum$mutation1[x]})
selection.write$gene2 <- sapply(1:nrow(selection.sum), function(x) if(selection.sum$direction[x] == 'above') {selection.sum$gene2[x]} else {selection.sum$gene1[x]})
selection.write$direction <- NULL
selection.write$q.value <- format(signif(selection.write$q.value,3), scientific = T)

selection.sig <- selection.write[which(selection.write$significant == 's'),]
gene.sig <- selection.sig[which(!(selection.sig$set %in% gene.sets$id)), c('set', 'mutation1', 'mutation2', 'disease', 'q.value')]
set.sig <- selection.sig[which(selection.sig$set %in% gene.sets$id), c('set', 'mutation1', 'gene1', 'mutation2', 'gene2', 'disease', 'q.value')]

write.table(gene.sig, file = paste0(wrk_dir, 'analysis/data_summary/Selection_Inference_Genes.csv'),
            sep = ',', row.names = F)
write.table(set.sig, file = paste0(wrk_dir, 'analysis/data_summary/Selection_Inference_Gene_Sets.csv'),
            sep = ',', row.names = F)

####################################################################################
# relative risks

pairwise.selection <- c()
for(gene in selection.genes){
  print(gene)
  if(gene %in% gene.sets$id){
    set.genes <- gene.sets[match(gene, gene.sets$id),-1][which(!is.na(gene.sets[match(gene, gene.sets$id),-1]))]
    mutations <- unname(unlist(sapply(set.genes, function(x) ListMutations(genes = x, threshold = freq.threshold))))
  } else {
    mutations <- ListMutations(gene, freq.threshold)
  }
  diseases <- unique(set.diseases$disease[which(set.diseases$set == gene)])
  for(disease in diseases){
    print(disease)
    res <- cbind(expand.grid(mutation1 = mutations, mutation2 = mutations), 
                 expand.grid(mutation1.id = 1:length(mutations), mutation2.id = 1:length(mutations), 
                             set = gene, 
                             disease = disease))
    ids <- metadata$patient.ID[which(metadata$disease.group == disease)]
    mut.sub <- mutation.aa.context.tab[match(mutations, rownames(mutation.aa.context.tab)),
                                       match(ids, colnames(mutation.aa.context.tab))]
    # exclude cases with more than one mutation from the set
    mut.sub <- mut.sub[, which(colSums(mut.sub) == 1)]
    res$gene1 <- hits.all$hugo.gene[match(res$mutation1, hits.all$aa.context)]
    res$gene2 <- hits.all$hugo.gene[match(res$mutation2, hits.all$aa.context)]
    res <- res[which(res$mutation1 != res$mutation2),]
    tmp <- unname(t(as.data.frame(lapply(1:nrow(res), function(x)
      RelativeRisk(mut1 = res$mutation1[x], 
                        mut2 = res$mutation2[x], 
                        mut.sub = mut.sub,
                        bootstrap.iterations = 100)))))
    res <- cbind(res, tmp)
    colnames(res)[(ncol(res) - 2):ncol(res)] <- c('rr', 'rr.lb', 'rr.ub')
    res$mutated.samples <- ncol(mut.sub)
    pairwise.selection <- rbind(pairwise.selection, res[,c('mutation1', 'mutation2', 'mutated.samples', 'set', 'disease', 'rr', 'rr.lb', 'rr.ub')])
  }
}

pairwise.selection$gene1 <- hits.all$hugo.gene[match(pairwise.selection$mutation1, hits.all$aa.context)]
pairwise.selection$gene2 <- hits.all$hugo.gene[match(pairwise.selection$mutation2, hits.all$aa.context)]

save(list = c('pairwise.selection'),
     file = paste(wrk_dir, 'R_processed_data/', 'PairwiseSelection.Rfile', sep = ''))


####################################################################################

# unitary analysis
load(file = paste(wrk_dir, 'R_processed_data/', 'UnitarySelection.Rfile', sep = ''))
unitary.selection <- c()
for(gene in selection.genes){
  print(gene)
  if(gene %in% gene.sets$id){
    set.genes <- gene.sets[match(gene, gene.sets$id),-1][which(!is.na(gene.sets[match(gene, gene.sets$id),-1]))]
    mutations <- unname(unlist(sapply(set.genes, function(x) ListMutations(genes = x, threshold = freq.threshold))))
  } else {
    mutations <- ListMutations(gene, freq.threshold)
  }
  diseases <- unique(set.diseases$disease[which(set.diseases$set == gene)])
  for(disease in diseases){
    print(disease)
    
    ids <- metadata$patient.ID[which(metadata$disease.group == disease)]
    mut.sub <- mutation.aa.context.tab[match(mutations, rownames(mutation.aa.context.tab)),
                                       match(ids, colnames(mutation.aa.context.tab))]
    
    res <- AnalyseSelection(mut.sub)
    res$set <- gene
    res$disease <- disease
    unitary.selection <- rbind(unitary.selection, res)
  }
}
save(list = c('unitary.selection'),
     file = paste(wrk_dir, 'R_processed_data/', 'UnitarySelection.Rfile', sep = ''))

################################################################################################################################
# variation explained

disease.totals <- table(metadata$disease.group)
unitary.selection$total <- disease.totals[match(unitary.selection$disease, names(disease.totals))]
unitary.selection$mutation.frequency <- unitary.selection$n / unitary.selection$total

res <- summary(lm(mutation.frequency ~ mean.p, data = unitary.selection))

explanation <- set.diseases
explanation$coefficient <- NA
explanation$r.squared <- NA
explanation$p.value <- NA
for(i in 1:nrow(explanation)){
  print(i)
  sub <- unitary.selection[which((unitary.selection$set == explanation$set[i]) & (unitary.selection$disease == explanation$disease[i])),]
  res <- summary(lm(n ~ mean.p, data = sub))
  explanation$coefficient[i] <- res$coefficients[2,1]
  explanation$p.value[i] <- res$coefficients[2,4]
  explanation$r.squared[i] <- res$r.squared
}
explanation <- explanation[order(explanation$p.value),]
explanation$stat <- sapply(1:nrow(explanation), function(x) (explanation$p.value[x] * nrow(explanation)) / length(which(explanation$p.value <= explanation$p.value[x])))
explanation$q.value <- sapply(1:nrow(explanation), function(x) min(explanation$stat[x:nrow(explanation)]))

sub <- explanation[(which(explanation$coefficient > 0)),]
(mean(sub$r.squared))
