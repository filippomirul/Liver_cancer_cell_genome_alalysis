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

setwd("C:\\Users\\filoa\\Desktop\\ResourceBioinfo\\R_Bioinfo_Res\\Project")
load("Liver_hepatocellular_carcinoma.RData")

#### Task 2 ####

ensembl <- useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl", mirror = "useast")

query <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name", "entrezgene_id"),
               filters = c("ensembl_gene_id"),
               values = list(c(r_anno_df$ensembl_gene_id)),
               mart = ensembl)

# Selecting only protein coding

query <- query[which(query$gene_biotype == "protein_coding"),]

r_anno_df <- r_anno_df[is.element(r_anno_df$ensembl_gene_id, query$ensembl_gene_id),]
raw_counts_df <- raw_counts_df[is.element(rownames(raw_counts_df), query$ensembl_gene_id),]

#### task 3 ####

count_thr <- 20
repl_thr <- 5

# creating a vector of rows that we need to remove
filter_vec <- apply(raw_counts_df, 1,
                    function(y) max(by(y, c_anno_df$condition, function(x) sum(x>=count_thr))))

# remove the rows from annotations and counts
filter_counts_df <- raw_counts_df[filter_vec>=repl_thr,]
filter_anno_df <- r_anno_df[rownames(filter_counts_df),]

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
DEGs <- as.data.frame(topTags(edge, n=20, p.value = 0.25, sort.by = "logFC"))
DEGs <- as.data.frame(topTags(edge, n=20000))
DEGs$class <- "="

# add labels for up (+) and down (-) regulations for certain parameters for -valued and FDR
DEGs$class[which(DEGs$logCPM > 1 & DEGs$logFC > 1.5 & DEGs$FDR < 0.25)] = "+"
DEGs$class[which(DEGs$logCPM > 1 & DEGs$logFC < -1.5 & DEGs$FDR < 0.25)] = "-"
DEGs <- DEGs[order(DEGs$logFC,decreasing = T),]

## Display an volcano plot
input_df <- DEGs
xlabel <- "log2 FC conrtol vs case"
ylabel <- "-log10 p-value"

par(fig=c(0,1,0,1), mar=c(4,4,1,2), mgp=c(2, 0.75, 0))	
plot(input_df$logFC, -log(input_df$PValue,base=10), xlab = xlabel, ylab = ylabel, 
     col = ifelse(input_df$class=="=","grey70", ifelse(input_df$class == "+", "indianred", "olivedrab4")), 
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
heatmap(as.matrix(CPM_table[temp, ]),
        ColSideColors = cols,cexCol = 0.5,margins = c(12,4),col=pal,cexRow = 0.2)
#### Task 4 ####

ensembl <- useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")
convert <- getBM(attributes=c("ensembl_gene_id","entrezgene_id","external_gene_name"),
                 filters=c("ensembl_gene_id"), 
                 values=DEGs$ensembl_gene_id,
                 mart = ensembl)

DEGs <- merge(DEGs, convert, by="ensembl_gene_id")
anyNA(DEGs)
# removing NAs
DEGs <- DEGs[which(!is.na(DEGs$entrezgene_id)),]
# removing duplicates
DEGs <- DEGs[-which(duplicated(DEGs$entrezgene_id)),]
# refining dataframe
DEGs <- DEGs %>% mutate(external_gene_name = external_gene_name.x) %>%
  select(-external_gene_name.y, -external_gene_name.x)

## Splitting up and down regulated ones

up_DEGs <- DEGs[DEGs$class == "+",]
down_DEGs <- DEGs[DEGs$class == "-",]

# Gene onthology analysis
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

barplot(up_ego_BP,showCategory=10)
barplot(down_ego_BP,showCategory=10)

dotplot(up_ego_BP, showCategory=20)
dotplot(down_ego_BP, showCategory = 20)

heatplot(up_ego_BP, showCategory = 6)
heatplot(down_ego_BP, showCategory = 6)

x2 <- pairwise_termsim(up_ego_BP) 
emapplot(x2)
x2 <- pairwise_termsim(down_ego_BP) 
emapplot(x2)

# MF
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

barplot(up_ego_MF,showCategory=10)
barplot(down_ego_MF,showCategory=10)

dotplot(up_ego_MF, showCategory=20)
dotplot(down_ego_MF, showCategory = 20)

heatplot(up_ego_MF, showCategory = 6)
heatplot(down_ego_MF, showCategory = 6)

x2 <- pairwise_termsim(up_ego_MF) 
emapplot(x2)
x2 <- pairwise_termsim(down_ego_MF) 
emapplot(x2)

# KEGG
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
         pathway.id = "hsa05202", 
         species = "human")

#### Task 6 ####

sequences <- getSequence(id = up_DEGs$ensembl_gene_id, mart = ensembl,
                         type = "ensembl_gene_id", upstream = 500, seqType = '5utr')
data("PWMLogn.hg19.MotifDb.Hsap")
seq <- lapply(sequences$cdna, function(x) DNAString(x))
enrichedTF <- motifEnrichment(seq[1:20], PWMLogn.hg19.MotifDb.Hsap, score = "affinity")
report <-  groupReport(enrichedTF)
report
plot(report[1:5])

#### Task 7 ####

#ID2
mdb.human <-  subset(MotifDb, organism =='Hsapiens' & geneSymbol == "ID2")
PWM <-  toPWM(as.list(mdb.human))
names(PWM) <-  sapply(names(PWM),function(x) strsplit(x,"-")[[1]][3])
scores <- motifScores(sequences = seq[1:5], PWM, raw.scores = T)

ecdf <-  motifEcdf(PWM, organism = "hg19", quick = TRUE)
threshold <-  log2(quantile(ecdf$ID2, 0.9975))
threshold
plotMotifScores(scores, sel.motifs="ID2", col = c("red","green","blue"), cutoff = threshold)

#### task 8 ####

tfs <- report$target[1:10]
scores <- data.frame(c(1:length(seq)))
#scores <- data.frame(c(1:5))


for (i in (1:5)){
  tfs_motif <- subset(MotifDb, organism == "Hsapiens" & geneSymbol == tfs[i])
  PWM <- toPWM(as.list(tfs_motif))
  names(PWM) <-  sapply(names(PWM),function(x) strsplit(x,"-")[[1]][3])

  ecdf <- motifEcdf(PWM, organism = "hg19", quick = T)
  threshold <-  log2(quantile(ecdf[[tfs[i]]], 0.9975))
  name <- tfs[i]
  
  sco <- motifScores(sequences = seq , PWM, cutoff = threshold)
  sco <-  as.data.frame(apply(sco, 1, sum))
  colnames(sco) <- name
  
  scores <- cbind(scores, sco)
  
  #fraction_9975 <- append(fraction_9975, length(which(apply(scores,1,sum) > 0)) /100)
}
 scores <- scores[,2:ncol(scores)]

result <- logical(nrow(scores))

# Loop through each row of the dataframe
for (i in 1:nrow(scores)) {
  # Check if any value in the row is equal to 0
  if (any(scores[i, ] < 3)) {
    result[i] <- FALSE
  } else {
    result[i] <- TRUE
  }
}



#### task 9 ####

up_names <- up_DEGs$external_gene_name
write.table(up_DEGs$external_gene_name,"up_DEGs.txt", sep="\n", row.names = F)
links <- read.delim("string_interactions_short.tsv")

#### Task 10 ####


ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
query <- getBM(attributes = c("external_gene_name","ensembl_gene_id","description",
                              "gene_biotype","start_position","end_position",
                              "chromosome_name","strand"),
               filters = c("ensembl_gene_id"), 
               values = up_DEGs[,1],
               mart = ensembl)

query <-  unique(query[,c(1,3:6)])

## Create the network
net <- igraph::graph_from_data_frame(d = links, vertices = query, directed = FALSE)







