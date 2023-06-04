#### Libraries ####
library(fgsea)
library(org.Hs.eg.db)
library(biomaRt)
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(DOSE)
library(pathview)
library(tidyverse)
library(dplyr)
library(edgeR)
library(limma)
library(GenomicFeatures)
library(ggplot2)
library(stringr)
library(tidyverse)
library(MotifDb)
library(seqLogo)
library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)
library(igraph)

#### Data ####
my_dir <- ""
setwd(my_dir)
load("Liver_hepatocellular_carcinoma.RData")
ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "hsapiens_gene_ensembl",
                      mirror = "useast")

#### Task 2 ####

query <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name", "entrezgene_id"),
               filters = c("ensembl_gene_id"),
               values = list(c(r_anno_df$ensembl_gene_id)),
               mart = ensembl)

# Selecting only protein coding

query <- query[which(query$gene_biotype == "protein_coding"),]

r_anno_df <- r_anno_df[is.element(r_anno_df$ensembl_gene_id,
                                  query$ensembl_gene_id),]
raw_counts_df <- raw_counts_df[is.element(rownames(raw_counts_df),
                                          query$ensembl_gene_id),]
raw_counts_df <- raw_counts_df + 1


#### task 3 ####

count_thr <- 20
repl_thr <- 5

# creating a vector of rows that we need to remove
filter_vec <- apply(raw_counts_df, 1,
                    function(y) max(by(y, c_anno_df$condition, function(x) sum(x >= count_thr))))

# remove the rows from annotations and counts
filter_counts_df <- raw_counts_df[filter_vec >= repl_thr,]
filter_anno_df <- r_anno_df[rownames(filter_counts_df), ]

edge_c <- DGEList(counts = filter_counts_df, group = c_anno_df$condition, 
                  samples = c_anno_df, genes = filter_anno_df) 

# normalization with Redge [TMM method]
edge_n <- calcNormFactors(edge_c, method = "TMM")
# create normalize expression values table
CPM_table <- as.data.frame(round(cpm(edge_n),2))

# design creation
design <- model.matrix(~ 0 + group, data = edge_n$samples)
colnames(design) <- levels(edge_n$samples$group)
rownames(design) <- edge_n$samples$sample

# calculate dispersion and fit with edgeR (necessary for differential expression analysis)
edge_d <- estimateDisp(edge_n, design)
edge_f <- glmQLFit(edge_d, design)

# definition of the contrast (conditions to be compared)
contrast <- makeContrasts("case-control", levels = design) 

# fit the model with generalized linear models
edge <- glmQLFTest(edge_f, contrast = contrast)
DEGs <- as.data.frame(topTags(edge, n = 20000, p.value = 0.01, sort.by = "logFC"))
DEGs$class <- "="

# add labels for up (+) and down (-) regulations for certain parameters for -valued and FDR
DEGs$class[which(DEGs$logCPM > 1 & DEGs$logFC > 1.5 )] = "+"
DEGs$class[which(DEGs$logCPM > 1 & DEGs$logFC < -1.5 )] = "-"
DEGs <- DEGs[order(DEGs$logFC, decreasing = T),]

input_df <- DEGs
xlabel <- "log2 FC conrtol vs case"
ylabel <- "-log10 p-value"

par(fig=c(0,1,0,1), mar=c(4,4,1,2), mgp=c(2, 0.75, 0))	
plot(input_df$logFC, -log(input_df$PValue,base=10), xlab = xlabel, ylab = ylabel, 
     col = ifelse(input_df$class=="=","grey85", ifelse(input_df$class == "+", "indianred", "olivedrab4")), 
     pch = 20, frame.plot = TRUE, cex = 0.8, main = "Volcano plot")
abline(v=0,lty=2,col="grey20")

# creation of vector of colors
cols <- c()
for (i in 1 : 100)
  if (c_anno_df$condition[i] == "control"){
    cols <- append(cols, "chartreuse4")
  } else {
    cols <- append(cols, "burlywood3")
  }

# select color palette
pal <- c("blue","white","red")
pal <- colorRampPalette(pal)(100)

# select disregulated genes
disregulated <- rownames(DEGs[DEGs$class != "=", ])

# create vector of binary values to select disregulated genes
temp <- rownames(CPM_table)%in%disregulated
#plot heatmap
heatmap(scale(log2(as.matrix(CPM_table[temp, ]))),
        ColSideColors = cols,cexCol = 0.5,margins = c(12,4),col=pal,cexRow = 0.2)
#### Task 4 ####

convert <- getBM(attributes=c("ensembl_gene_id","entrezgene_id"),
                 filters=c("ensembl_gene_id"), 
                 values=DEGs$ensembl_gene_id,
                 mart = ensembl)

DEGs <- merge(DEGs, convert, by = "ensembl_gene_id")
anyNA(DEGs)
# removing NAs
DEGs <- DEGs[which(!is.na(DEGs$entrezgene_id)),]
# removing duplicates
DEGs <- DEGs[-which(duplicated(DEGs$entrezgene_id)),]

## Splitting up and down regulated ones

up_DEGs <- DEGs[DEGs$class == "+",]
down_DEGs <- DEGs[DEGs$class == "-",]

# Gene ontology analysis
# BP
up_ego_BP <- enrichGO(gene = up_DEGs$external_gene_name,
                      OrgDb = org.Hs.eg.db,
                      keyType = 'SYMBOL',
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

down_ego_BP <- enrichGO(gene = down_DEGs$external_gene_name,
                        OrgDb = org.Hs.eg.db,
                        keyType = 'SYMBOL',
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
dotplot(up_ego_BP, showCategory = 10)
dotplot(down_ego_BP, showCategory = 10)

heatplot(up_ego_BP, showCategory = 6)
heatplot(down_ego_BP, showCategory = 6)

x2 <- pairwise_termsim(up_ego_BP) 
emapplot(x2)
x2 <- pairwise_termsim(down_ego_BP) 
emapplot(x2)

up_ego_MF <- enrichGO(gene = up_DEGs$external_gene_name,
                      OrgDb = org.Hs.eg.db,
                      keyType = 'SYMBOL',
                      ont = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

down_ego_MF <- enrichGO(gene = down_DEGs$external_gene_name,
                        OrgDb = org.Hs.eg.db,
                        keyType = 'SYMBOL',
                        ont = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

dotplot(up_ego_MF, showCategory = 10)
dotplot(down_ego_MF, showCategory = 10)

heatplot(up_ego_MF, showCategory = 6)
heatplot(down_ego_MF, showCategory = 6)

x2 <- pairwise_termsim(up_ego_MF) 
emapplot(x2)
x2 <- pairwise_termsim(down_ego_MF) 
emapplot(x2)

up_eKEGG <- enrichKEGG(gene = up_DEGs$entrezgene_id,
                       organism = 'human',
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1)

down_eKEGG <- enrichKEGG(gene = down_DEGs$entrezgene_id,
                         organism = 'human',
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.1)
#### Task 5 ####

FC <- up_DEGs$logFC
names(FC) <- up_DEGs$entrezgene_id
pathview(gene.data = FC, 
         pathway.id = "hsa04110",  # hsa04512 (ECM receptor interaction)
         species = "human")
#### Task 6 ####

sequences <- getSequence(id = up_DEGs$ensembl_gene_id, mart = ensembl,
                         type = "ensembl_gene_id", upstream = 500, seqType = "gene_flank")

seq <- lapply(sequences[,1], function(x) DNAString(x))

#data("PWMLogn.hg19.MotifDb.Hsap")
# enrichedTF <- motifEnrichment(seq, PWMLogn.hg19.MotifDb.Hsap, score = "affinity")
# save(enrichedTF, file = "enriched_transcription.RData")
load("enriched_transcription.RData")

report <- groupReport(enrichedTF)
# we have selected only those above the 0.975 percentage of the data
report <- report[report$p.value < 0.025,]
report[1:5, ]

#### Task 7 ####

mdb.human <-  subset(MotifDb, organism =='Hsapiens' & geneSymbol == "PGAM2")
PWM <-  toPWM(as.list(mdb.human))
names(PWM) <-  sapply(names(PWM),function(x) strsplit(x,"-")[[1]][3])
scores <- motifScores(sequences = seq, PWM, raw.scores = T)

ecdf <-  motifEcdf(PWM, organism = "hg19", quick = TRUE)
threshold <-  log2(quantile(ecdf$PGAM2, 0.9975))

plotMotifScores(scores[1:8], sel.motifs = "PGAM2",
                col = c("red","green","blue"),
                cutoff = threshold)

#### task 8 ####

tfs <- report$target
# Some TFs gave us problem during the analysis (all the ones with the terms UW.M),
#  because they cannot be found in the MotifDb, so we had to cut them out.
tfs <- tfs[!grepl("^UW.M", tfs)]
scores <- data.frame(c(1:length(seq)))


# for (i in (1:length(tfs))){
# 
#   tfs_motif <- subset(MotifDb, organism == "Hsapiens" & geneSymbol == tfs[i])
#   PWM <- toPWM(as.list(tfs_motif))
#   names(PWM) <-  sapply(names(PWM),function(x) strsplit(x,"-")[[1]][3])
# 
#   ecdf <- motifEcdf(PWM, organism = "hg19", quick = T)
#   threshold <-  log2(quantile(ecdf[[tfs[i]]], 0.9975))
#   name <- tfs[i]
# 
#   sco <- motifScores(sequences = seq , PWM, cutoff = threshold)
#   sco <-  as.data.frame(apply(sco, 1, sum))
#   colnames(sco) <- name
# 
#   scores <- cbind(scores, sco)
# }
# save(scores, file = "Scores_table.RData")


load("Scores_table.RData")
# The first column is the index 
scores <- as.matrix(scores[, 2:length(scores)])
head(scores[, 1:10])


# filter out TFs that do not have any score greater than 0 for any of the sequences
columns <- c()
for (col in (1:ncol(scores))){
  if (sum(scores[,col]) < 1){
    columns <- append(columns, col)
  }
}
# removing the zero columns
if (length(columns) > 0){
  scores <- scores[ ,-columns]  
}

# initialized a treshold for the cut, in this case the median, but can be increased (has been used median instead of mean because more stable).
# dist is a sub set of all the scores to confront the line (sequence) distrbution of TFs
threshold <- 0.5 # 50%
dist <- c()

for (i in(1:ncol(scores))){
  point <- sample(scores[,i], 150, replace = F)
  dist <- append(dist, point)
}

q <- quantile(dist, probs = c(threshold))
# the loop return a boolean to filter each line (gene up-strem sequence),
# if the distribution of the row has a higher value of the threshold the row is kept.

result <- logical(nrow(scores))
for (i in 1:nrow(scores)) {
  # 
  if (quantile(scores[i,], probs = c(threshold)) < q) {
    result[i] <- FALSE
  } else {
    result[i] <- TRUE
  }
}
summary(result)

#### task 9 ####

top_enrich <- sequences$ensembl_gene_id[result]

up_names <- up_DEGs$external_gene_name[which(top_enrich%in%up_DEGs$ensembl_gene_id)]

## the command below produce a txt file with the gene of interest that have been used on STRING to retrive the TSV file.

lapply(up_names, write,my_dir , append=TRUE, ncolumns=1000)

links <- read.delim("string_interactions_short.tsv")

#### Task 10 ####

query <- getBM(attributes = c("external_gene_name","ensembl_gene_id","description",
                              "gene_biotype","start_position","end_position",
                              "chromosome_name","strand"),
               filters = c("ensembl_gene_id"), 
               values = top_enrich,
               mart = ensembl)

query <-  unique(query[,c(1,3:6)]) %>% arrange(external_gene_name)
index <- links$node2%in%query$external_gene_name
links <- links[index,]
index <- links$X.node1%in%query$external_gene_name
links <- links[index,]

## Create the network
net <- graph_from_data_frame(d = links, vertices = query, directed = FALSE)
net_clean <- simplify(net)

# Plot graph
plot(net_clean)

c <- components(net_clean, mode =  c("weak","strong"))
## find max and position
max(c$csize)
net_work <- which.max(c$csize)
net.c <-  induced_subgraph(net_clean,V(net_clean)[which(c$membership == net_work)])
plot(net.c,
     edge.width = 5,
     vertex.color = "orange",
     vertex.size = 10,
     vertex.frame.color = "darkgray",
     vertex.label.color = "black", 
     vertex.label.cex = 0.7,
     edge.curved = 0.1)

ceb <- cluster_edge_betweenness(net.c) 
plot(ceb, net.c, vertex.size = 10)

# retrieving the names from each gene list
hclust_obj <- as.hclust(ceb)
clusters <- cutree(hclust_obj, k = 11)
cluster_names <- lapply(unique(clusters), function(cluster) {
  names(clusters)[clusters == cluster]
})
for (i in seq_along(cluster_names)) {
  cat("Cluster", i, ": ", paste(cluster_names[[i]], collapse = ", "), "\n")
}
