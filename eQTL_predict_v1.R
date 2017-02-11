source("predict_functions.R")
source("predict_expression_funcs.R")
require(GOSim)
require(randomForest)

gm.combined.mat <- as.matrix(read.table("./input/CRM_matrix.txt", sep="\t", header=TRUE))
tf.names <- colnames(gm.combined.mat)[1:(ncol(gm.combined.mat) -3)]

# load eQTL data
dimas.proxy <- as.matrix(read.table("./input/eQTL_data.txt", sep="\t", header=TRUE))
LCL.proxy <- dimas.proxy[which(dimas.proxy[,"Tissue"]=="L"),]

# load gene expression data
haem.express <- as.matrix(read.table("./input/expression_data.txt", sep="\t", header=TRUE))
all.samples <- 2:ncol(haem.express)

# gene annotation data 
gene.pos <- as.matrix(read.table("./input/gene_positions_hg18.txt", header=TRUE, sep="\t"))


# load FAIRE data
faire.gene <- as.matrix(read.table("./input/GM12878_FAIREseq_hg18.pk", sep="\t", header=TRUE))

# load CTCF data
gm.ctcf <-  as.matrix(read.table("./input/GM12878_CTCF.bed", sep="\t", header=TRUE))


# unique identifier for each combination of TF binding sites
gm.combined.mat <- getPattern(gm.combined.mat, ncol(gm.combined.mat) -3)

# find regions with co-localization of CRM and eQTL SNPs
crm.snps <- overlap.eqtl(gm.combined.mat, LCL.proxy)

# measure co-expression between TFs and genes
co.express <- regress.expression(crm.snps, haem.express, all.samples, gene.pos=gene.pos)

# number of crms to analyze
samps <- nrow(crm.snps)


gene.crms <- genes.to.crms(gm.combined.mat, gene.pos)
faire.gene <- genes.to.crms(faire.gene, gene.pos)

gene.sim <- GO.similarity(crm.snps, gene.pos, tf.names)

insulators <- map.ctcf(crm.snps, gene.pos, gm.ctcf)

score.mat <- cooccur.score(gm.combined.mat[,1:length(tf.names)])

save.image("GM12878_29TFs_Dimas_prediction.RData")

features.mat <- feature.select(gene.preds=co.express, crm.mat=crm.snps, samps=samps[-which(is.na(co.express))], faire.gene=faire.gene, gene.crms=gene.crms, gene.sim=gene.sim, insulators=insulators, score.mat=score.mat)

# combines 

