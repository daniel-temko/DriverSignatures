
TestChangeWilcox <- function(presence.table, signature.table, metadata, mutation, group, process, group.by, alternate){
  presence.table.sub <- presence.table[ , which(metadata[[group.by]] == group)]
  signature.scores.table.sub <- signature.table[ , which(metadata[[group.by]] == group)]
  in.indices <- which(presence.table.sub[match(mutation, rownames(presence.table.sub)), ] > 0)
  
  if(length(in.indices) == ncol(presence.table.sub)){
    return(c(NA, NA))
  }
  #process.sub <- signature.tissue.df[which(signature.tissue.df$short.name == short.name), ]
  sig <- wilcox.test(as.numeric(signature.scores.table.sub[match(process, rownames(signature.scores.table.sub)),in.indices]),
                     as.numeric(signature.scores.table.sub[match(process, rownames(signature.scores.table.sub)),-in.indices]),
                     alternative = alternate,
                     conf.int = T)
  effect.size <- sig$estimate
  p.value <- sig$p.value
  return(c(effect.size, p.value))
}

RankSignatureEffectOnContext <- function(signature, context, disease){
  sigs <- signature.tissue.df$variable[which(signature.tissue.df$signature.group == disease)]
  sigs <- sigs[which(sigs != 'Other')]
  rank <- match(signature, sigs[order(prc.proven.even[match(context, processes),match(sigs, colnames(prc.proven.even))], decreasing = T)])
  return(rank)
}

ContextPropInDisease <- function(signature, context, disease){
  sigs <- signature.tissue.df$variable[which(signature.tissue.df$signature.group == disease)]
  sigs <- sigs[which((sigs != 'Other') & (sigs != signature))]
  av.level <- mean(as.numeric(prc.proven.even[match(context, processes), match(sigs, colnames(prc.proven.even))]))
  return(av.level)
}


####################################################################################################################################
# driver - signature associations

## Create the data-frame to use for tests
n.base <- 3
data <- AggregateOccurencesByCancer(mutation.aa.context.tab, metadata, 'disease.group')
data$mutation <- rownames(data)
hits <- melt(data, id.vars = 'mutation')
hits <- hits[which(hits$value >= n.base), ]
colnames(hits)[match('variable', colnames(hits))] <- 'signature.group'
hits <- merge(hits, signature.tissue.df, by = 'signature.group')
hits$signature.group <- droplevels(hits$signature.group)
colnames(hits)[match('variable', colnames(hits))] <- 'signature'
hits <- hits[which(hits$signature != 'Other'),]
hits$context <- hits.all$context[match(hits$mutation, hits.all$aa.context)]
hits$context.prop.in.sig <- sapply(1:nrow(hits), function(x) prc.proven.even[match(hits$context[x], processes),match(hits$signature[x], colnames(prc.proven.even))])
hits$context.prop.in.disease <- sapply(1:nrow(hits), function(x) ContextPropInDisease(signature = as.character(hits$signature[x]), context = hits$context[x], disease = as.character(hits$signature.group[x])))
hits$context.lfc <- log((hits$context.prop.in.sig + (1/96)) / (hits$context.prop.in.disease + (1/96)))
hits <- hits[which(hits$context.lfc > 0),]

# test significance
load(file = paste(wrk_dir, 'R_processed_data/', 'AARes.Rfile', sep = ''))
aa.res <- t(as.data.frame(lapply(1:nrow(hits), function(x) TestChangeWilcox(presence.table = mutation.aa.context.tab,
                                                                            signature.table = signature.scores.norm, 
                                                                            metadata = metadata, 
                                                                            mutation = as.character(hits$mutation[x]), 
                                                                            group = as.character(hits$signature.group[x]), 
                                                                            process = hits$signature[x], 
                                                                            group.by = 'disease.group',
                                                                            alternate = 'greater'))))
save(list = c('aa.res'),
     file = paste(wrk_dir, 'R_processed_data/', 'AARes.Rfile', sep = ''))
hits <- cbind(hits, aa.res)
colnames(hits)[((ncol(hits) - 1):ncol(hits))] <- c('effect.size', 'p.value')

# analyse number of hits compared to mutation filter
n.vals <- n.base:20
n.hits <- c()
for(i in n.base:20){
  n <- i
  aa.hits <- hits[which(hits$value >= n),]
  aa.hits <- aa.hits[order(aa.hits$p.value),]
  aa.hits$stat <- sapply(1:nrow(aa.hits), function(x) (aa.hits$p.value[x] * nrow(aa.hits)) / length(which(aa.hits$p.value <= aa.hits$p.value[x])))
  aa.hits$q.value <- sapply(1:nrow(aa.hits), function(x) min(aa.hits$stat[x:nrow(aa.hits)]))
  fdr <- 0.05
  n.hits <- c(n.hits, length(which((aa.hits$q.value < fdr))))
}
pdf(file = paste(wrk_dirm 'analysis/data_summary/Hits.pdf', sep = ''))
plot(n.vals, n.hits)
dev.off()
n <- 4
aa.hits <- hits[which(hits$value >= n),]
aa.hits <- aa.hits[order(aa.hits$p.value),]
aa.hits$stat <- sapply(1:nrow(aa.hits), function(x) (aa.hits$p.value[x] * nrow(aa.hits)) / length(which(aa.hits$p.value <= aa.hits$p.value[x])))
aa.hits$q.value <- sapply(1:nrow(aa.hits), function(x) min(aa.hits$stat[x:nrow(aa.hits)]))

####################################################################################################################################
# Summarise the results for downstream analysis

aa.hits$mutation.frequency <- aa.hits$value / sapply(1:nrow(aa.hits), function(x) length(which(metadata$disease.group == aa.hits$signature.group[x])))
aa.hits$gene <- hits.all$hugo.gene[match(aa.hits$mutation, hits.all$aa.context)]
aa.hits$classification <- hits.all$classification[match(aa.hits$mutation, hits.all$aa.context)]
aa.hits$annotation <- signature.annotations$association[match(aa.hits$signature, signature.annotations$process)]
aa.hits$signature.type <- sapply(1:nrow(aa.hits), function(x) signature.annotations$type[match(aa.hits$signature[x], signature.annotations$process)])
aa.hits$significant <- sapply(aa.hits$q.value, function(x) if (x < fdr) {'s'} else {'ns'})
aa.hits$direction <- sapply(aa.hits$effect.size, function(x) if(x > 0){ 'pos'} else {'neg'})

aa.sum <- aa.hits[,c('mutation', 'gene', 'signature', 'annotation', 'signature.type', 'signature.group', 'classification', 'value', 'mutation.frequency', 'context', 'effect.size', 'context.prop.in.sig', 'context.prop.in.disease', 'context.lfc', 'q.value', 'significant', 'direction')]
save(list = c('aa.sum'),
     file = paste(wrk_dir, 'R_processed_data/', 'AASum.Rfile', sep = ''))


aa.sig <- aa.sum[which(aa.sum$significant == 's'), ]
aa.sig$samples.in.tumour.type <- sapply(1:nrow(aa.sig), function(x) length(which(metadata$disease.group == aa.sig$signature.group[x])))
aa.write <- aa.sig[order(aa.sig$signature), c('mutation', 'gene', 'signature', 'annotation', 'signature.group', 'context', 'value', 'samples.in.tumour.type', 'mutation.frequency', 'q.value')]
aa.write$context <- sapply(1:nrow(aa.write), function(x) TransformContext(aa.write$context[x]))
aa.write$mutation.frequency <- sprintf("%.2f", aa.write$mutation.frequency)
aa.write$q.value <- format(signif(aa.write$q.value,3), scientific = T)
aa.write$mutation <- hits.all$short.name[match(aa.write$mutation, hits.all$aa.context)]
colnames(aa.write)[match(c('signature.group', 'value', 'context', 'mutation.frequency'), colnames(aa.write))] <- c('disease', 'mutation.count', 'causal.channel', 'mutation.frequency.in.tumour.type')

# write aa results lits
write.table(aa.write, file = paste0(wrk_dir, 'analysis/data_summary/Driver_Signature_Associations.csv'), sep = ',', row.names = F)
