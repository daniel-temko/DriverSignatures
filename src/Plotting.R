require(ggrepel)
require(RColorBrewer)
require(grid)
require(gridExtra)

#width = 6.9
#height = 4.9
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector[c(7, 13)] <- col_vector[26 + c(3, 4)]

ContextPropInDisease <- function(signature, context, disease){
  sigs <- signature.tissue.df$variable[which(signature.tissue.df$signature.group == disease)]
  sigs <- sigs[which((sigs != 'Other') & (sigs != signature))]
  av.level <- mean(as.numeric(prc.proven.even[match(context, processes), match(sigs, colnames(prc.proven.even))]))
  return(av.level)
}

out_dir <- ""
data_dir <- ""

################################################################################################################################
# Figure 1 - Overview figure

width = 6.6
height = 4.9
cols <- rep(c('blue', 'black', 'red', 'grey', 'green', 'pink'), each = 16)

# toy signature
sig.s <- prc.proven.exome[,15] + runif(96)*0.01
sig.s[c(43, 71)] <- sig.s[c(71, 43)]
sig.s <- sig.s / sum(sig.s)
heights1 <- sig.s
sig.av <- sapply(1:nrow(prc.proven.exome), function(x) mean(as.numeric(cbind(prc.proven.exome[x,c(5,11:12)], sig.s[x]))))
heights2 <- sig.av
y.max <- max(c(heights1, heights2)) * 1.1
cols1 <- sapply(cols, function(x) adjustcolor(x, alpha.f = 0.2))
cols2 <- cols
context <- 'TC C.G'
context.lab <- 'CTG>CCG'
process.id <- match(context, processes)

pdf(file = paste(out_dir, '/Figure.1B.pdf', sep = ''),
    width = width, height = height)
par(mar = c(4,6,4,2))
bar <- barplot(heights1, 
               col = cols1, ylim = c(0, y.max), las = 2, main = 'Mutational Signatures', cex.axis = 1.5,
               cex.main = 2, border = F, yaxt = 'n')
par(new = T)
bar <- barplot(heights2, col = cols2, ylim = c(0, y.max), yaxt = 'n', 
               cex.main = 1.5, border = F)
legend('topleft', legend = c('Signature S', 'Average Signature'), pch = 15, 
       col = c(adjustcolor('black', alpha.f = 0.5), 'black'), cex = 1.5)
axis(side = 1, at = bar[process.id], labels = context.lab, cex.axis = 2)
axis(side = 2, at = seq(0, y.max, 0.05), labels = F)
title(ylab = 'Probability', cex.lab = 2)
dev.off()

# toy boxplot
mut.df <- data.frame(type = 'mut', signature = rnorm(100, 3))
wt.df <- data.frame(type = 'wt', signature = rnorm(100, 2))
sig.sub <- rbind(mut.df, wt.df)
sig.sub$signature <- sig.sub$signature / max(sig.sub$signature)
sig.sub$type <- factor(sig.sub$type, levels = c('wt', 'mut'))

pdf(file = paste(out_dir, '/Figure.1C.pdf', sep = ''),
    width = width, height = height)
par(mar = c(4,8,5,2))
box <- boxplot(signature ~ type, data = sig.sub, col = c('grey', 'green'), pch = 16, ylim = c(0, 1), 
               xaxt = 'n', yaxt = 'n', ylab = 'Signature S', cex.lab = 2,
               main = 'Signature S in Samples \nwith Mutation A',
               cex.main = 2)
axis(side = 1, at = c(1, 2),
     labels = c('wt', 'Mutation A'), padj = 1, cex.axis = 2)
axis(side = 2, at = seq(0, 1, 0.2), labels = F)
dev.off()

# second toy signature
heights1 <- sig.av
heights1[39] <- heights1[71] * (2/5)
#y.max <- max(heights1) * 1.1
contexts <- c('CT G.G', 'TC C.G')
context.labs <- c('GCG>GTG', 'CTG>CCG')
process.ids <- match(contexts, processes)

pdf(file = paste(out_dir, '/Figure.1E.pdf', sep = ''),
    width = width, height = height)
par(mar = c(4,6,4,2))
bar <- barplot(heights1, 
               col = cols, ylim = c(0, y.max), las = 2, main = 'Mutational Signatures', cex.axis = 1.5,
               cex.main = 2, border = F, yaxt = 'n')
legend('topleft', legend = c('Sample Signature'), pch = 15, 
       col = c('black'), cex = 1.5)
axis(side = 1, at = bar[process.ids[1]], labels = context.labs[1], cex.axis = 2)
axis(side = 1, at = bar[process.ids[2]], labels = F, cex.axis = 2)
axis(side = 2, at = seq(0, y.max, 0.05), labels = F)
title(ylab = 'Probability', cex.lab = 2)
dev.off()


width = 6.6
height = 5.9

## Frequency barplot and relative risk
# frequency barplot
prob.a <- 0.05
freq.a <- 0.055
data.a <- matrix(c(prob.a, freq.a - prob.a, freq.a, 0), nrow = 2)
prob.b <- 0.02
freq.b <- 0.08
data.b <- matrix(c(prob.b, freq.b - prob.b, freq.b, 0), nrow = 2)

#data <- matrix(c(prob.a, freq.a, prob.b, freq.b), nrow = 2)
y.max <- max(c(data.a, data.b)) * 1.2
col.a <- 'green'
col.b <- 'red'
pdf(file = paste(out_dir, '/Figure.1F.pdf', sep = ''),
    width = width, height = height)
par(mar = c(10,7,5,2))
x <- barplot(cbind(data.a, matrix(rep(0,8), nrow = 2)), 
             ylim = c(0, y.max), space = 0, 
             col = c(col.a, adjustcolor(col.a, alpha.f = 0.2)), axes = F,
             names.arg = c('Expected', 'Observed', '', '', 'Expected', 'Observed'), main = 'Expected and Observed Mutations',
             cex.main = 2, cex.names = 1.5, las = 2)
par(new = T)
barplot(cbind(matrix(rep(0,8), nrow = 2), data.b), 
        ylim = c(0, y.max), space = 0,
        col= c(col.b, adjustcolor(col.b, alpha.f = 0.2)), axes = F)
axis(side = 1, at = c(sum(x[1:2])/2, sum(x[5:6])/2), 
     labels = c('Mutation A', 'Mutation B'), tick = F, line = 7,
     cex.axis = 1.5)
axis(side = 2, labels = F)
#axis(side = 2)
#axis(side = 4)
#title(ylab = 'Probability')
dev.off()

# relative risk
risk.a <- freq.a - prob.a
risk.b <- freq.b - prob.b
pdf(file = paste(out_dir, '/Figure.1G.pdf', sep = ''),
    width = width, height = height)
par(mar = c(10,7,5,2))
x <- barplot(c(0, risk.a, rep(0,2), risk.b, 0), ylim = c(0, y.max), ylab = 'Risk',
             col = c(0, adjustcolor(col.a, alpha.f = 0.2), rep(0, 2), adjustcolor(col.b, alpha.f = 0.2), 0),
             names.arg = '', main = 'Mutation Risk', 
             cex.main = 2, cex.lab = 2,
             cex.names = 1.5, yaxt = 'n')
lines(x = c(0, max(x)), y = c(0, 0))
axis(side = 1, at = x[c(2,5)], 
     labels = c('Mutation A', 'Mutation B'), tick = F,
     cex.axis = 2)
axis(side = 2, labels = F)
#names.arg = c('Mutation A', rep('', 3), 'Mutation B', rep('', 3), 'Mutation C')
dev.off()

################################################################################################################################
# Figure 2 - Association figure

load(file = paste(data_dir, 'AAPower.Rfile', sep = ''))
aa.sig <- aa.power[which(aa.power$significant == 's'), ]

width = 6.9
height = 4.9

## correlation plots
figure.nums <- paste('2', LETTERS[1:6], sep = '')
sigs <- c('Signature.10', 'Signature.4', 'Signature.26', 'Signature.2', 'Signature.1', 'Signature.15')
sig.annotations <- c('POLE', 'Smoking', 'MMR Defects', 'APOBEC', 'Ageing', 'MMR Defects')
diseases <- c('Uterine Carcinoma', 'Lung Adeno', 'Uterine Carcinoma', 'Breast', 'Colorectum', 'Stomach')
disease.names <- c('UC', 'LUAD', 'UC', 'BRCA', 'CRC', 'Stomach')
muts <- c('PTEN R130Q (CT T.G)', 'KRAS G12C (CA C.A)', 'KRAS G12D (CT A.C)', 'PIK3CA E545K (CT T.A)', 'APC R213X (CT A.G)', 'FBXW7 R465C (CT G.G)')
context.names <- paste('(', c('TCG>TTG', 'CCA>CAA' ,'ACC>ATC', 'TCA>TTA', 'ACG>ATG', 'GCG>GTG'), ')', sep = '')

sig.names <- sapply(sigs, function(x) paste(strsplit(x, '\\.')[[1]], collapse = ' '))
gene.names <- hits.all$hugo.gene[match(muts, hits.all$aa.context)]
aa.change.names <- hits.all$residue[match(muts, hits.all$aa.context)]


for(i in 1:length(figure.nums)){
  disease.ids <- metadata$patient.ID[which(metadata$disease.group == diseases[i])]
  mut.ids <- colnames(mutation.aa.context.tab)[which(mutation.aa.context.tab[match(muts[i], rownames(mutation.aa.context.tab)),] > 0)]
  sig.sub <- data.frame(signature = signature.scores.norm[match(sigs[i], rownames(signature.scores.norm)), match(disease.ids, colnames(signature.scores.norm))],
                        type = sapply(disease.ids, function(x) if(x %in% mut.ids) { 'mut' } else { 'wt' }))
  sig.sub$type <- factor(sig.sub$type, levels = c('wt', 'mut'))
  q.val <- aa.sig$q.value[which((aa.sig$signature.group == diseases[i]) & (aa.sig$mutation == muts[i]) & (aa.sig$signature == sigs[i]))]
  
  pdf(file = paste(out_dir, '/Figure.', figure.nums[i], '.pdf', sep = ''),
      width = width, height = height)
  par(mar = c(5,5,5,2))
  box <- boxplot(signature ~ type, data = sig.sub, pch = 16, col = c('blue', 'red'), ylim = c(0, 1.4), 
                 xaxt = 'n', yaxt = 'n', ylab = sig.names[i], cex.lab = 1.5,
                 main = paste(paste(gene.names[i], aa.change.names[i], context.names[i]), ' and\n', sig.names[i], ' in ', disease.names[i], sep = ''), 
                 cex.main = 1.5)
  axis(side = 1, at = c(1, 2), 
       labels = c('wt', paste(gene.names[i], aa.change.names[i], '\n', context.names[i])),
       padj = 1, cex.axis = 1.5)
  axis(side = 2, at = seq(0, 1, 0.2))
  
  # add significance
  highest.prop <- max(sig.sub$signature)
  text(x = c(1,2), y = rep(highest.prop + 0.32, 2), paste('N = ', table(sig.sub$type)), cex = 1.5)
  lines(x = c(1, 2), y = c(highest.prop + 0.1, highest.prop + 0.1))
  text(x = 1.5, y = highest.prop + 0.2, paste('Q=', format(signif(q.val, 2), scientific = T), sep = ''), cex = 1.5)
  dev.off()
}

################################################################################################################################
# Figure 3 - within-gene selection

load(file = paste(data_dir, 'PairwiseSelection.Rfile', sep = ''))

width = 6.9
height = 4.9

figure.diseases <- c('Pancreas', 'Breast', 'Glioma Low Grade', 'Lung Adeno')
figure.genes <- c('KRAS', 'PIK3CA', 'IDH1', 'KRAS')

color.map <- data.frame(gene = unique(figure.genes))
color.map$col <- col_vector[1:nrow(color.map)]

## barplots
figure.nums <- c('3A', '3B', '3C', '3D')
for(i in 1:length(figure.diseases)){
  gene <- figure.genes[i]
  disease <- figure.diseases[i]
  figure.num <- figure.nums[i]
  
  sub <- pairwise.selection[which((pairwise.selection$set == gene) & (pairwise.selection$disease == disease)),]
  sub$rr.length <- sub$rr.ub - sub$rr.lb
  
  mutations <- ListMutations(gene, freq.threshold)
  summary <- data.frame(mut = mutations)
  summary$confidence <- sapply(1:nrow(summary), function(x) 1 / mean(sub$rr.length[which(sub$mutation2 == as.character(summary$mut[x]))]))
  #summary$selection.evidence <- sapply(1:nrow(summary), function(x) length(which((sub$mutation1 == as.character(summary$mut[x])) & (sub$significant == 's'))))
  summary$short.name <- hits.all$short.name[match(summary$mut, hits.all$aa.context)]
  summary$residue <- hits.all$residue[match(summary$mut, hits.all$aa.context)]
  summary$context.name <- hits.all$context.name[match(summary$mut, hits.all$aa.context)]
  summary$label <- sapply(1:nrow(summary), function(x) if(length(which(summary$residue == summary$residue[x])) == 1) {summary$residue[x]} else {paste(summary$residue[x], ' (', summary$context[x], ')', sep = '')})
  ref <- as.character(summary$mut[(which.max(summary$confidence))])
  ref.name <- as.character(summary$label[match(ref, summary$mut)])
  
  sub <- pairwise.selection[which((pairwise.selection$set == gene) &
                               (pairwise.selection$disease == disease) &
                               (pairwise.selection$mutation2 == ref)),]
  sub <- rbind(sub, data.frame(mutation1 = ref, mutation2 = ref, mutated.samples = sub$mutated.samples[1], set = gene, disease = disease, rr = 1, rr.lb = 1, rr.ub = 1, gene1 = gene, gene2 = gene))
  sub <- sub[order(sub$rr),]
  sub$label <- summary$label[match(sub$mutation1, summary$mut)]
  sub$length.ci <- sub$rr.ub - sub$rr.lb
  
  y.max <- max(max(sub$rr.ub) * 1.1, 1.4)
  col <- color.map$col[match(gene, color.map$gene)]
  pdf(file = paste(out_dir, '/Figure.', figure.num, '.pdf', sep = ''),
      width = width, height = height)
  par(mar = c(13, 6, 4, 6))
  x <- barplot(sub$rr, ylim = c(0, y.max), names.arg = sub$label, las = 2, col = col,
               ylab = paste('RR vs', ref.name), main = paste(gene, 'in', disease), cex.lab = 1.5, 
               cex.names = 1.5, cex.main = 1.5)
  abline(h = 1, lty = 2, col = 'grey')
  text(median(x), 0.8 * y.max, labels = paste('N=', sub$mutated.samples[1], sep =''), cex = 1.5)
  for(j in 1:nrow(sub)){
    if(sub$length.ci[j] == 0){next}
    arrows(x[j,], sub$rr.lb[j], x[j,], sub$rr.ub[j], angle = 90, code = 3, length = 0.1)
  }
  dev.off()
}


################################################################################################################################
# Figure 4 - KRAS, BRAF and NRAS selection

width = 8
height = 5.6

load(file = paste(data_dir, 'PairwiseSelection.Rfile', sep = ''))

figure.sets <- rep('gs.1', 6)
figure.diseases <- c('Pancreas', 'Uterine Carcinoma', 'Thyroid', 'Melanoma', 'Colorectum', 'Lung Adeno')

# barplots
figure.nums <- c('4A', '4B', '4C', '4D', '4E', '4F')

for(i in 1:length(figure.diseases)){
  gene <- figure.sets[i]
  disease <- figure.diseases[i]
  figure.num <- figure.nums[i]
  
  sub <- pairwise.selection[which((pairwise.selection$set == gene) & (pairwise.selection$disease == disease)),]
  sub$rr.length <- sub$rr.ub - sub$rr.lb
  
  set.genes <- gene.sets[match(gene, gene.sets$id),-1][which(!is.na(gene.sets[match(gene, gene.sets$id),-1]))]
  mutations <- unname(unlist(sapply(set.genes, function(x) ListMutations(genes = x, threshold = freq.threshold))))
  
  summary <- data.frame(mut = mutations)
  summary$confidence <- sapply(1:nrow(summary), function(x) 1 / mean(sub$rr.length[which(sub$mutation2 == as.character(summary$mut[x]))]))
  #summary$selection.evidence <- sapply(1:nrow(summary), function(x) length(which((sub$mutation1 == as.character(summary$mut[x])) & (sub$significant == 's'))))
  summary$short.name <- hits.all$short.name[match(summary$mut, hits.all$aa.context)]
  summary$residue <- hits.all$residue[match(summary$mut, hits.all$aa.context)]
  summary$context.name <- hits.all$context.name[match(summary$mut, hits.all$aa.context)]
  summary$gene <- hits.all$hugo.gene[match(summary$mut, hits.all$aa.context)]
  summary$label <- sapply(1:nrow(summary), function(x) if(length(which((summary$residue == summary$residue[x]) & (summary$gene == summary$gene[x]))) == 1) {paste(summary$gene[x], summary$residue[x])} else {paste(summary$gene[x], ' ', summary$residue[x], ' (', summary$context[x], ')', sep = '')})
  ref <- as.character(summary$mut[which.max(summary$confidence)])
  ref.name <- as.character(summary$label[which.max(summary$confidence)])
  
  sub <- pairwise.selection[which((pairwise.selection$set == gene) &
                               (pairwise.selection$disease == disease) &
                               (pairwise.selection$mutation2 == ref)),]
  sub <- rbind(sub, data.frame(mutation1 = ref, mutation2 = ref, mutated.samples = sub$mutated.samples[1], set = gene, disease = disease, rr = 1, rr.lb = 1, rr.ub = 1, gene1 = gene, gene2 = gene))
  sub <- sub[order(sub$rr),]
  sub$gene <- summary$gene[match(sub$mutation1, summary$mut)]
  sub$label <- summary$label[match(sub$mutation1, summary$mut)]
  sub$length.ci <- sub$rr.ub - sub$rr.lb
  
  color.map <- data.frame(gene = unique(sort(sub$gene)))
  color.map$col <- col_vector[1:nrow(color.map)]
  sub$col <- color.map$col[match(sub$gene, color.map$gene)]
  
  max.y <- max(max(sub$rr.ub) * 1.1, 1.4)
  pdf(file = paste(out_dir, '/Figure.', figure.num, '.pdf', sep = ''),
      width = width, height = height)
  par(mar = c(13, 6, 6, 6))
  x <- barplot(sub$rr, names.arg = sub$label, ylim = c(0, max.y),
               col = as.character(sub$col), las = 2, main = disease, 
               ylab = paste('RR vs', ref.name), cex.lab = 1.5, cex.main = 1.5)
  text(median(x), 0.8 * max.y, labels = paste('N=', sub$mutated.samples[1], sep = ''), cex = 1.5)
  abline(h = 1, lty = 2, col = 'grey')
  if(i == 1){
    legend('topleft', pch = 15, legend = color.map$gene, col = color.map$col,
           cex = 1.5)
  }
  for(j in 1:nrow(sub)){
    if(sub$length.ci[j] == 0){next}
    arrows(x[j,], sub$rr.lb[j], x[j,], sub$rr.ub[j], angle = 90, code = 3, length = 0.1)
  }
  dev.off()
}

################################################################################################################################
# Figure 5 - APC and CTNNB1

width = 14
height = 5.6

load(file = paste(data_dir, 'PairwiseSelection.Rfile', sep = ''))

figure.sets <- rep('gs.2', 3)
figure.diseases <- c('Liver', 'Uterine Carcinoma', 'Colorectum')

# scatters
figure.nums <- c('5.A', '5.B', '5.C')

for(i in 1:length(figure.diseases)){
  gene <- figure.sets[i]
  disease <- figure.diseases[i]
  figure.num <- figure.nums[i]
  
  sub <- pairwise.selection[which((pairwise.selection$set == gene) & (pairwise.selection$disease == disease)),]
  sub$rr.length <- sub$rr.ub - sub$rr.lb
  
  set.genes <- gene.sets[match(gene, gene.sets$id),-1][which(!is.na(gene.sets[match(gene, gene.sets$id),-1]))]
  mutations <- unname(unlist(sapply(set.genes, function(x) ListMutations(genes = x, threshold = freq.threshold))))
  
  summary <- data.frame(mut = mutations)
  summary$confidence <- sapply(1:nrow(summary), function(x) 1 / mean(sub$rr.length[which(sub$mutation2 == as.character(summary$mut[x]))]))
  #summary$selection.evidence <- sapply(1:nrow(summary), function(x) length(which((sub$mutation1 == as.character(summary$mut[x])) & (sub$significant == 's'))))
  summary$short.name <- hits.all$short.name[match(summary$mut, hits.all$aa.context)]
  summary$residue <- hits.all$residue[match(summary$mut, hits.all$aa.context)]
  summary$context.name <- hits.all$context.name[match(summary$mut, hits.all$aa.context)]
  summary$gene <- hits.all$hugo.gene[match(summary$mut, hits.all$aa.context)]
  summary$label <- sapply(1:nrow(summary), function(x) if(length(which((summary$residue == summary$residue[x]) & (summary$gene == summary$gene[x]))) == 1) {paste(summary$gene[x], summary$residue[x])} else {paste(summary$gene[x], summary$residue[x], ' (', summary$context[x], ')', sep = '')})
  ref <- as.character(summary$mut[which.max(summary$confidence)])
  ref.name <- as.character(summary$label[which.max(summary$confidence)])
  
  sub <- pairwise.selection[which((pairwise.selection$set == gene) &
                               (pairwise.selection$disease == disease) &
                               (pairwise.selection$mutation2 == ref)),]
  sub <- rbind(sub, data.frame(mutation1 = ref, mutation2 = ref, mutated.samples = sub$mutated.samples[1], set = gene, disease = disease, rr = 1, rr.lb = 1, rr.ub = 1, gene1 = gene, gene2 = gene))
  sub <- sub[order(sub$rr),]
  sub$gene <- summary$gene[match(sub$mutation1, summary$mut)]
  sub$label <- summary$label[match(sub$mutation1, summary$mut)]
  sub$length.ci <- sub$rr.ub - sub$rr.lb
  
  color.map <- data.frame(gene = unique(sort(sub$gene)))
  color.map$col <- col_vector[1:nrow(color.map)]
  sub$col <- color.map$col[match(sub$gene, color.map$gene)]
  
  max.y <- max(max(sub$rr.ub) * 1.1, 1.2)
  pdf(file = paste(out_dir, '/Figure.', figure.num, '.pdf', sep = ''),
      width = width, height = height)
  par(mar = c(13, 6, 6, 6))
  x <- barplot(sub$rr, names.arg = sub$label, ylim = c(0, max.y),
               col = as.character(sub$col), las = 2, main = disease, 
               ylab = paste('RR vs', ref.name), cex.lab = 1.5, cex.main = 1.5)
  text(median(x), 0.8 * max.y, labels = paste('N=', sub$mutated.samples[1], sep = ''), cex = 1.5)
  abline(h = 1, lty = 2, col = 'grey')
  if(i == 1){
    legend('topleft', pch = 15, legend = color.map$gene, col = color.map$col,
           cex = 1.5)
  }
  for(j in 1:nrow(sub)){
    if(sub$length.ci[j] == 0){next}
    arrows(x[j,], sub$rr.lb[j], x[j,], sub$rr.ub[j], angle = 90, code = 3, length = 0.1)
  }
  dev.off()
}

########################################################################################
## Supplementary figures
########################################################################################

# Figures S1-3

load(file = paste(out_dir, 'Summary.Rfile', sep = ''))

#raster plots
width =  8.27 * 2
height = 11.69 * 2 
figure.nums <- 1:2

k = 0
for(j in 1:2){
  figure.name <- paste('Figure.S', figure.nums[j], sep = '')
  #figure.num <- paste('Figure.S1.', tolower(as.roman(j)), sep = '')
  #figure.name <- gsub('\\.S', ' S', figure.num)
  file = paste(out_dir, figure.name, '.pdf', sep = '')
  pdf(file = file, width = width, height = height)
  par(mfrow = c(5,3), mar = c(4,4,4,2), oma = c(10,10,10,10))
  for(i in 1:15){
    k = 15*(j - 1) + i
    plot(unique(scheme), summary[k,], type = 'l', lty = 1, xlim = c(0,96), ylim = c(0,1), xlab = '', ylab = '', 
         main = paste('Signature', k), cex.main = 2, cex.axis = 2, lwd = 2)
  }
  #title(main = figure.name, outer = T, cex.main = size.1)
  mtext('Mutations in Sample', side = 1, outer = T, cex = 2, line = 3)
  mtext('Classification Accuracy', side = 2, outer = T, cex = 2, line = 3)
  dev.off()
}

av.accuracy <- apply(summary, 2, mean)
figure.num <- 'Figure.S3'
figure.name <- gsub('\\.S', ' S', figure.num)
file = paste(out_dir, figure.num, '.pdf', sep = '')
pdf(file = file, width = 11.69, height = 8.27)
par(mar = c(4,4,4,2), oma = c(4,10,8,8))
plot(unique(scheme), av.accuracy, type = 'l', lty = 1, xlim = c(0,96), ylim = c(0,1), xlab = 'Mutations in Sample', ylab = 'Classification Accuracy', 
     main = 'Average Across Signatures', cex.main = 1.5, cex.lab = 1.2)
abline(v = 20, h = 0.8)
#title(main = figure.name, outer = T, cex.main = size.1)
dev.off()

########################################################################################

# Figure S4
load(file = paste(data_dir, 'AASum.Unbiased.Rfile', sep = ''))

(tab <- table(aa.sum$mechanistic.basis, aa.sum$significant))

test <- fisher.test(matrix(as.numeric(table(aa.sum$mechanistic.basis, aa.sum$significant)), nrow = 2))

pval <- format(signif(test$p.value,2), scientific = T)
#### CHECK HARDCODED FORMULAE AGAINST PVAL #####
#a <- expression(italic('P') * '=5.5e-5')
#### CHECK HARDCODED FORMULAE AGAINST PVAL #####

n.sig <- aa.sum$context.lfc[which(aa.sum$significant == 'ns')]
sig <- aa.sum$context.lfc[which(aa.sum$significant == 's')]
cols <- data.frame(col = c('blue', 'red', 'grey'), alpha = c(0.5, 0.5, 0.5))
cols$hex <- sapply(1:nrow(cols), function(x) do.call(rgb,as.list(c(col2rgb(cols$col[x])/256, cols$alpha[x]))))
pdf(file = paste(out_dir, '/Figure.S4.pdf', sep = ''),
    width = 7, height = 10)
par(mar = c(6,6,8,4))
breaks <- seq(-5,5,0.25)
hist(n.sig, col = cols$hex[3], breaks = breaks, 
     xlab = 'Log fold change causal channel', main = 'All tests', 
     cex.lab = 2, cex.axis = 2, cex.main = 2, lty = 'blank')
hist(sig, col = cols$hex[2] ,breaks = breaks, add = T, lty = 'blank')
abline(v = 0, lty = 2)
text(x = 2.8, 
     y = 120, labels = paste('P=', pval, sep = ''), cex = 1.7)
legend('topleft', pch = 15, legend = c('Association', 'No Association'), 
       col = c(cols$hex[2], cols$hex[3]), cex = 1.5)
dev.off()

########################################################################################

# Figures S5-15

width = 8.27 * 2
height = 11.69 * 2

load(file = paste(data_dir, 'PairwiseSelection.Rfile', sep = ''))

#figure.nums <- paste('S', 3:11, sep = '')
figure.genes <- selection.genes[1:9]
col <- col_vector[1]
rows <- 4
text.size <- 2
x.lab.size <- 2
y.lab.size <- 2
y.axis.size <- 1.5
main.size <- 2
x.mar <- c(8, 16, rep(8, 7))

figure.num <- 5
for(i in 1:length(figure.genes)){
  gene <- figure.genes[i]
  print(gene)
  mutations <- unique(pairwise.selection$mutation1[which(pairwise.selection$set == gene)])
  if(length(mutations) < 15){cols <- 2} else {cols <- 1}
  #figure.num <- figure.nums[i]
  
  diseases <- unique(set.diseases$disease[which(set.diseases$set == gene)])
  for(j in 1:length(diseases)){
    disease <- as.character(diseases[j])
    if((j %% (cols * rows)) == 1){
      figure.name <- paste('S', figure.num, sep = '')
      pdf(file = paste(out_dir, '/Figure.', figure.name, '.pdf', sep = ''),
          width = width, height = height)
      par(mfrow = c(rows,cols), oma = c(2,4,6,4), mar = c(x.mar[i], 6, 4, 6))
      figure.num <- figure.num + 1
    }
    print(disease)
    sub <- pairwise.selection[which((pairwise.selection$set == gene) & (pairwise.selection$disease == disease)),]
    sub$rr.length <- sub$rr.ub - sub$rr.lb
    
    summary <- data.frame(mut = mutations)
    #summary$mean.rr <- sapply(1:nrow(summary), function(x) mean(sub$rr[which(sub$mutation2 == as.character(summary$mut[x]))], na.rm = T))
    summary$confidence <- sapply(1:nrow(summary), function(x) 1 / mean(sub$rr.length[which(sub$mutation2 == as.character(summary$mut[x]))]))
    summary$short.name <- hits.all$short.name[match(summary$mut, hits.all$aa.context)]
    summary$residue <- hits.all$residue[match(summary$mut, hits.all$aa.context)]
    summary$context.name <- hits.all$context.name[match(summary$mut, hits.all$aa.context)]
    summary$label <- sapply(1:nrow(summary), function(x) if(length(which(summary$residue == summary$residue[x])) == 1) {summary$residue[x]} else {paste(summary$residue[x], ' (', summary$context[x], ')', sep = '')})
    ref <- as.character(summary$mut[which.max(summary$confidence)])
    ref.name <- as.character(summary$label[which.max(summary$confidence)])
    
    sub <- pairwise.selection[which((pairwise.selection$set == gene) &
                                 (pairwise.selection$disease == disease) &
                                 (pairwise.selection$mutation2 == ref)),]
    sub <- rbind(sub, data.frame(mutation1 = ref, mutation2 = ref, mutated.samples = sub$mutated.samples[1], set = gene, disease = disease, rr = 1, rr.lb = 1, rr.ub = 1, gene1 = gene, gene2 = gene))
    sub <- sub[order(sub$rr),]
    sub$label <- summary$label[match(sub$mutation1, summary$mut)]
    sub$length.ci <- sub$rr.ub - sub$rr.lb
    
    finite.ubs <- sub$rr.ub[which(!is.infinite(sub$rr.ub))]
    y.max <- max(max(sub$rr) * 1.1, max(finite.ubs) * 1.1, 1.2)
    x <- barplot(sub$rr, ylim = c(0, y.max), col = col, names.arg = sub$label, 
                 las = 2, ylab = paste('RR vs', ref.name), main = paste(gene, 'in', disease), cex.lab = y.lab.size, cex.names = x.lab.size, 
                 cex.main = main.size, cex.axis = y.axis.size)
    text(median(x), 0.8 * y.max, labels = paste('N=', sub$mutated.samples[1], sep = ''), cex = text.size)
    abline(h = 1, lty = 2, col = 'grey')
    for(k in 1:nrow(sub)){
      if(sub$length.ci[k] == 0){next}
      arrows(x[k,], sub$rr.lb[k], x[k,], min(sub$rr.ub[k], y.max * 1.2), angle = 90, code = 3, length = 0.1)
    }
    if(((j %% (cols*rows)) == 0) | (j == length(diseases))) {
      #title(main = paste('Figure', figure.name), outer = T)
      dev.off()
    }
  }
}

########################################################################################
# Figures S16-S20

width = 8.27 * 3
height = 11.69 * 3
rows <- 5
text.size <- 3
x.lab.size <- 2
y.lab.size <- 3
y.axis.size <- 2
main.size <- 3
x.mar <- rep(23, 3)
#x.mar <- c(8, 16, rep(8, 7))


figure.genes <- selection.genes[-(1:9)]

figure.num <- 16
for(i in 1:length(figure.genes)){
  gene <- figure.genes[i]
  print(gene)
  
  set.name <- gene.sets$name[match(gene, gene.sets$id)]
  mutations <- unique(pairwise.selection$mutation1[which(pairwise.selection$set == gene)])
  if(length(mutations) < 15){cols <- 2} else {cols <- 1}
  
  diseases <- unique(set.diseases$disease[which(set.diseases$set == gene)])
  panel.num <- 1
  for(j in 1:length(diseases)){
    disease <- as.character(diseases[j])
    if((j %% (cols * rows)) == 1){
      figure.name <- paste('S', figure.num, sep = '')
      pdf(file = paste(out_dir, '/Figure.', figure.name, '.pdf', sep = ''),
          width = width, height = height)
      par(mfrow = c(rows,cols), oma = c(2,4,6,4), mar = c(x.mar[i], 6, 4, 6))
      figure.num <- figure.num + 1
      print(disease)
    }
    
    sub <- pairwise.selection[which((pairwise.selection$set == gene) & (pairwise.selection$disease == disease)),]
    sub$rr.length <- sub$rr.ub - sub$rr.lb
    
    summary <- data.frame(mut = mutations)
    #summary$mean.rr <- sapply(1:nrow(summary), function(x) mean(sub$rr[which(sub$mutation2 == as.character(summary$mut[x]))], na.rm = T))
    summary$confidence <- sapply(1:nrow(summary), function(x) 1 / mean(sub$rr.length[which(sub$mutation2 == as.character(summary$mut[x]))]))
    summary$short.name <- hits.all$short.name[match(summary$mut, hits.all$aa.context)]
    summary$residue <- hits.all$residue[match(summary$mut, hits.all$aa.context)]
    summary$context.name <- hits.all$context.name[match(summary$mut, hits.all$aa.context)]
    summary$gene <- hits.all$hugo.gene[match(summary$mut, hits.all$aa.context)]
    summary$label <- sapply(1:nrow(summary), function(x) if(length(which((summary$residue == summary$residue[x]) & (summary$gene == summary$gene[x]))) == 1) {paste(summary$gene[x], summary$residue[x])} else {paste(summary$gene[x], summary$residue[x], ' (', summary$context[x], ')', sep = '')})
    ref <- as.character(summary$mut[which.max(summary$confidence)])
    ref.name <- as.character(summary$label[which.max(summary$confidence)])
    
    sub <- pairwise.selection[which((pairwise.selection$set == gene) &
                                 (pairwise.selection$disease == disease) &
                                 (pairwise.selection$mutation2 == ref)),]
    sub <- rbind(sub, data.frame(mutation1 = ref, mutation2 = ref, mutated.samples = sub$mutated.samples[1], set = gene, disease = disease, rr = 1, rr.lb = 1, rr.ub = 1, gene1 = gene, gene2 = gene))
    sub <- sub[order(sub$rr),]
    sub$gene <- summary$gene[match(sub$mutation1, summary$mut)]
    sub$label <- summary$label[match(sub$mutation1, summary$mut)]
    sub$length.ci <- sub$rr.ub - sub$rr.lb
    
    color.map <- data.frame(gene = unique(sort(sub$gene)))
    color.map$col <- col_vector[1:nrow(color.map)]
    sub$col <- color.map$col[match(sub$gene, color.map$gene)]
    
    finite.ubs <- sub$rr.ub[which(!is.infinite(sub$rr.ub))]
    y.max <- max(max(sub$rr) * 1.1, max(finite.ubs) * 1.1, 1.2)
    x <- barplot(sub$rr, ylim = c(0, y.max), names.arg = sub$label, las = 2, col = as.character(sub$col),
                 ylab = paste('RR vs', ref.name), main = paste(set.name, 'in', disease), cex.lab = y.lab.size, cex.names = x.lab.size, 
                 cex.main = main.size, cex.axis = y.axis.size, yaxt = 'n')
    text(median(x), 0.8 * y.max, labels = paste('N=', sub$mutated.samples[1], sep = ''), cex = text.size)
    legend('topleft', legend = color.map$gene, col = color.map$col, pch = 15, cex = 3)
    axis(side = 2, at = seq(0, y.max, 0.5), cex.axis = y.axis.size)
    abline(h = 1, lty = 2, col = 'grey')
    for(k in 1:nrow(sub)){
      if(sub$length.ci[k] == 0){next}
      arrows(x[k,], sub$rr.lb[k], x[k,], min(sub$rr.ub[k], y.max * 1.2), angle = 90, code = 3, length = 0.1)
    }
    #title(main = paste('Figure', figure.name), outer = T)
    if(((j %% (cols*rows)) == 0) | (j == length(diseases))) {
      #title(main = paste('Figure', figure.name), outer = T)
      dev.off()
    }
  }
}

########################################################################################
# Figures S21-S22

width = 6.6
height = 4.9

load(file = paste(data_dir, 'SelectionSum.Rfile', sep = ''))
selection.write <- selection.sum
selection.write$mutation1 <- sapply(1:nrow(selection.sum), function(x) if(selection.sum$direction[x] == 'above') {selection.sum$mutation1[x]} else {selection.sum$mutation2[x]})
selection.write$gene1 <- sapply(1:nrow(selection.sum), function(x) if(selection.sum$direction[x] == 'above') {selection.sum$gene1[x]} else {selection.sum$gene2[x]})
selection.write$mutation2 <- sapply(1:nrow(selection.sum), function(x) if(selection.sum$direction[x] == 'above') {selection.sum$mutation2[x]} else {selection.sum$mutation1[x]})
selection.write$gene2 <- sapply(1:nrow(selection.sum), function(x) if(selection.sum$direction[x] == 'above') {selection.sum$gene2[x]} else {selection.sum$gene1[x]})
selection.write$direction <- NULL
selection.write$mutation.disease <- paste(selection.write$mutation1, selection.write$disease)

set.write.sum <- selection.write[which(selection.write$set %in% gene.sets$id),]
set.write.sig <- set.write.sum[which(set.write.sum$significant == 's'),]

set <- 'gs.1'
sub <- set.write.sig[which(set.write.sig$set == set),]
sub <- sub[match(unique(sub$mutation.disease), sub$mutation.disease),]
tab <- table(sub$gene1, droplevels(sub$disease))
tmp <- melt(tab, varnames = c('Gene', 'Disease'))
figure.num <- 'S21'
pdf(file = paste(out_dir, '/Figure.', figure.num, '.pdf', sep = ''),
    width = width, height = height, onefile = F)
par(mfrow = c(2,1))
f <- ggplot(tmp, aes(Disease, Gene)) + 
  geom_tile(aes(fill = value), colour = 'black') + 
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0, na.value = 'white') +
  geom_text(aes(label = value)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.text = element_text(size = 16)) + 
  theme(axis.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) + 
  theme(legend.title = element_text(size = 16))
print(f)
dev.off()


set <- 'gs.2'
sub <- set.write.sig[which(set.write.sig$set == set),]
sub <- sub[match(unique(sub$mutation.disease), sub$mutation.disease),]
tab <- table(sub$gene1, droplevels(sub$disease))
tmp <- melt(tab, varnames = c('Gene', 'Disease'))
figure.num <- 'S22'
pdf(file = paste(out_dir, '/Figure.', figure.num, '.pdf', sep = ''),
    width = width, height = height, onefile = F)
f <- ggplot(tmp, aes(Disease, Gene)) + 
  geom_tile(aes(fill = value), colour = 'black') + 
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0, na.value = 'white') +
  geom_text(aes(label = value)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.text = element_text(size = 16)) + 
  theme(axis.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) + 
  theme(legend.title = element_text(size = 16))
print(f)
dev.off()


########################################################################################
# Figures S23-S33

load(file = paste(data_dir, 'UnitarySelection.Rfile', sep = ''))

figure.genes <- selection.genes[1:9]
figure.nums <- paste('S', 17:25, sep = '')

width = 8.27 * 2
height = 11.69 * 2
rows <- 4
cols <- 2

figure.num <- 23
for(i in 1:length(figure.genes)){
  gene <- figure.genes[i]
  print(gene)
  #figure.num <- figure.nums[i]
  diseases <- as.character(set.diseases$disease[which(set.diseases$set == gene)])
  for(j in 1:((ceil(length(diseases) / (rows * cols))) *  (rows * cols))){
    panel.num <- j %% (rows * cols)
    if(panel.num == 1){
      figure.name <- paste('Figure.S', figure.num, sep = '')
      pdf(file = paste(out_dir, figure.name, '.pdf', sep = ''),
          height = height, width = width)
      par(oma = c(6,6,6,6))
      figure.num <- figure.num + 1
    }
    if(j > length(diseases)){
      f <- rectGrob(gp=gpar(col="white"))
    } else {
      disease <- diseases[j]
      sub <- unitary.selection[which((unitary.selection$set == gene) & 
                                     (unitary.selection$disease == disease)),]
      sub$short.name <- hits.all$short.name[match(sub$mut, hits.all$aa.context)]
      sub$residue <- hits.all$residue[match(sub$mut, hits.all$aa.context)]
      sub$context.name <- hits.all$context.name[match(sub$mut, hits.all$aa.context)]
      sub$label <- sapply(1:nrow(sub), function(x) if(length(which(sub$residue == sub$residue[x])) == 1) {sub$residue[x]} else {paste(sub$residue[x], ' (', sub$context[x], ')', sep = '')})
      labels <- sub[which(sub$n > 0),]
      res <- summary(lm(n ~ mean.p, data = sub))
      
      x.max <- max(sub$mean.p) * 1.2
      y.max <- max(sub$n) * 1.2
      f <- ggplot(sub, aes(mean.p, n)) +
        coord_cartesian(xlim = c(0, x.max), ylim = c(0, y.max)) +
        geom_point(aes(color = gene), size = 2) +
        geom_abline(intercept = res$coefficients[1,1], slope = res$coefficients[2,1], linetype = 2, color = 'grey', lwd = 0.5) + 
        ggtitle(paste(gene, ' in ', disease, sep = '')) +
        xlab('Mean Probability') +
        ylab('Frequency') +
        geom_text_repel(data = labels, aes(label = label), force = 1, segment.alpha = 0.5, size = 5) +
        theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
        theme(axis.text=element_text(size=18),
               axis.title=element_text(size=20,face="bold")) + 
        theme(plot.title = element_text(size = 20, face = "bold")) +
        theme(legend.text=element_text(size= 20), legend.title=element_blank())
    }
    assign(paste('g', panel.num, sep = ''), f)
    if(panel.num == 0){
      grid.arrange(g1, g2, g3, g4, g5, g6, g7, g0, ncol = 2)
      dev.off()
    }
  }
}

########################################################################################
# Figures S34-S38

load(file = paste(data_dir, 'UnitarySelection.Rfile', sep = ''))

figure.genes <- selection.genes[-(1:9)]
figure.nums <- paste('S', 26:28, sep = '')

width = 8.27 * 2
height = 11.69 * 2
rows <- 5
cols <- 1

figure.num <- 34
for(i in 1:length(figure.nums)){
  gene <- figure.genes[i]
  gene.name <- gene.sets$name[match(gene, gene.sets$id)]
  print(gene.name)
  diseases <- as.character(set.diseases$disease[which(set.diseases$set == gene)])
  for(j in 1:((ceil(length(diseases) /  (rows * cols))) * (rows * cols))){
    panel.num <- j %% (rows * cols)
    if(panel.num == 1){
      figure.name <- paste('Figure.S', figure.num, sep = '')
      pdf(file = paste(out_dir, figure.name, '.pdf', sep = ''),
          height = height, width = width)
      par(oma = c(6,6,6,6))
      figure.num <- figure.num + 1
    }
    if(j > length(diseases)){
      f <- rectGrob(gp=gpar(col="white"))
    } else {
      disease <- diseases[j]
      sub <- unitary.selection[which((unitary.selection$set == gene) & 
                                       (unitary.selection$disease == disease)),]
      sub$grad <- sub$n / sub$mean.p
      sub$short.name <- hits.all$short.name[match(sub$mut, hits.all$aa.context)]
      sub$residue <- hits.all$residue[match(sub$mut, hits.all$aa.context)]
      sub$context.name <- hits.all$context.name[match(sub$mut, hits.all$aa.context)]
      sub$gene <- hits.all$hugo.gene[match(sub$mut, hits.all$aa.context)]
      sub$label <- sapply(1:nrow(sub), function(x) if(length(which((sub$residue == sub$residue[x]) & (sub$gene == sub$gene[x])) ) == 1) {sub$residue[x]} else {paste(sub$residue[x], ' (', sub$context[x], ')', sep = '')})
      labels <- sub[which(match(sub$mut, sub$mut[order(sub$grad, decreasing = T)]) <= 4),]
      res <- summary(lm(n ~ mean.p, data = sub))
      
      x.max <- max(sub$mean.p) * 1.2
      y.max <- max(sub$n) * 1.2
      f <- ggplot(sub, aes(mean.p, n)) +
        coord_cartesian(xlim = c(0, x.max), ylim = c(0, y.max)) +
        geom_point(aes(color = gene), size = 2) +
        geom_abline(intercept = res$coefficients[1,1], slope = res$coefficients[2,1], linetype = 2, color = 'grey', lwd = 0.5) + 
        ggtitle(paste(gene.name, ' in ', disease, sep = '')) +
        xlab('Mean Probability') +
        ylab('Frequency') +
        geom_text_repel(data = labels, aes(label = label), force = 1, segment.alpha = 1, size = 5) +
        theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
        theme(axis.text=element_text(size=18),
              axis.title=element_text(size=20,face="bold")) + 
        theme(plot.title = element_text(size = 20, face = "bold")) +
        theme(legend.text=element_text(size= 18), legend.title = element_blank())
    }
    assign(paste('g', panel.num, sep = ''), f)
    if(panel.num == 0){
      grid.arrange(g1, g2, g3, g4, g0, ncol = 2)
      dev.off()
    }
  }
}
