

Analyze of the first experiment: NCms10058
==========================================

Cluster and annotate in the shell (not in R)
--------------------------------------------


```bash
LIBRARY=NCms10058_1
BAMFILES=../Moirai/NCms10058_1.CAGEscan_short-reads.20150625154711/properly_paired_rmdup/*bam

level1.py -o $LIBRARY.l1.gz -f 66 -F 516 $BAMFILES 
level2.py -t 0 -o $LIBRARY.l2.gz $LIBRARY.l1.gz

function osc2bed {
  zcat $1 |
    grep -v \# |
    sed 1d |
    awk '{OFS="\t"}{print $2, $3, $4, "l1", "1000", $5}'
}

function bed2annot {
  bedtools intersect -a $1 -b ../annotation/annot.bed -s -loj |
    awk '{OFS="\t"}{print $1":"$2"-"$3$6,$10}' | 
    bedtools groupby -g 1 -c 2 -o collapse
}

function bed2symbols {
  bedtools intersect -a $1 -b ../annotation/gencode.v14.annotation.genes.bed -s -loj |
    awk '{OFS="\t"}{print $1":"$2"-"$3$6,$10}' | 
    bedtools groupby -g 1 -c 2 -o distinct
}

osc2bed $LIBRARY.l2.gz | tee $LIBRARY.l2.bed | bed2annot - > $LIBRARY.l2.annot
bed2symbols $LIBRARY.l2.bed > $LIBRARY.l2.genes
```

```
## Opening NCms10058_1.l1.gz
```

Analysis with R
---------------

### Configuration


```r
library(oscR)        #  See https://github.com/charles-plessy/oscR for oscR.
library(smallCAGEqc) # See https://github.com/charles-plessy/smallCAGEqc for smallCAGEqc.
library(vegan)
```

```
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.0-10
```

```r
library(ggplot2)

stopifnot(
    packageVersion("oscR") >= "0.1.1"
  , packageVersion("smallCAGEqc") > "0.10.0"
)

LIBRARY <- "NCms10058_1"
```

### Load data


```r
l2_NCki <- read.osc(paste(LIBRARY,'l2','gz',sep='.'), drop.coord=T, drop.norm=T)
colnames(l2_NCki) <- sub('raw.NCms10058_1.','NCki_',colnames(l2_NCki))
colSums(l2_NCki)  
```

```
##    NCki_HeLa_PS_A    NCki_HeLa_PS_B    NCki_HeLa_PS_C NCki_HeLa_RanN6_A NCki_HeLa_RanN6_B NCki_HeLa_RanN6_C 
##             11800             13969             22764             14137             13556             10430 
##    NCki_THP1_PS_A    NCki_THP1_PS_B    NCki_THP1_PS_C NCki_THP1_RanN6_A NCki_THP1_RanN6_B NCki_THP1_RanN6_C 
##             15157             15453             13092              8708             14536             17122
```

### Normalization number of read per sample : l2.sub ; libs$genes.sub

In all the 3 libraries used, one contain only few reads tags. The smallest one has 8,708 counts. In order to make meaningful comparisons, all of them are subsapled to 8700 counts.


```r
set.seed(1)
l2.sub1 <- t(rrarefy(t(l2_NCki),min(8700)))
colSums(l2.sub1)
```

```
##    NCki_HeLa_PS_A    NCki_HeLa_PS_B    NCki_HeLa_PS_C NCki_HeLa_RanN6_A NCki_HeLa_RanN6_B NCki_HeLa_RanN6_C 
##              8700              8700              8700              8700              8700              8700 
##    NCki_THP1_PS_A    NCki_THP1_PS_B    NCki_THP1_PS_C NCki_THP1_RanN6_A NCki_THP1_RanN6_B NCki_THP1_RanN6_C 
##              8700              8700              8700              8700              8700              8700
```

### Moirai statistics

Load the QC data produced by the Moirai workflow with which the libraries were processed. Sort in the same way as the l1 and l2 tables, to allow for easy addition of columns.


```r
libs <- loadMoiraiStats(multiplex = "NCms10058_1.multiplex.txt", summary = "../Moirai/NCms10058_1.CAGEscan_short-reads.20150625154711/text/summary.txt", pipeline = "CAGEscan_short-reads")
```

### Number of clusters

Count the number of unique L2 clusters per libraries after subsampling, and add
this to the QC table.  Each subsampling will give a different result, but the
mean result can be calculated by using the `rarefy` function at the same scale
as the subsampling.


```r
libs["l2.sub1"]     <- colSums(l2.sub1 > 0)
libs["l2.sub1.exp"] <- rarefy(t(l2_NCki), min(colSums(l2_NCki)))
```

### Richness

Richness should also be calculated on the whole data.


```r
libs["r100.l2"] <- rarefy(t(l2_NCki),100)
boxplot(data=libs, r100.l2 ~ group, ylim=c(80,100), las=1)
```

![](commandes_clean_files/figure-html/NCki_boxplot-1.png) 

### Hierarchical annotation

Differences of sampling will not bias distort the distribution of reads between annotations, so the non-subsampled library is used here.


```r
annot.l2 <- read.table(paste(LIBRARY,'l2','annot',sep='.'), head=F, col.names=c('id', 'feature'), row.names=1)
annot.l2 <- hierarchAnnot(annot.l2)

rownames(libs) <- sub("HeLa", "NCki_HeLa", rownames(libs))
rownames(libs) <- sub("THP1", "NCki_THP1", rownames(libs))

libs <- cbind(libs, t(rowsum(l2_NCki,  annot.l2[,'class']))) 
libs$samplename <- sub('HeLa', 'NCki_HeLa', libs$samplename)
libs$samplename <- sub('THP1', 'NCki_THP1', libs$samplename)
```

### Gene symbols used normalisation data


```r
genesymbols <- read.table(paste(LIBRARY,'l2','genes',sep='.'), col.names=c("cluster","symbol"), stringsAsFactors=FALSE)
rownames(genesymbols) <- genesymbols$cluster

g2 <- rowsum(l2_NCki, genesymbols$symbol)
countSymbols <- countSymbols(g2)

libs[colnames(l2_NCki),"genes"] <- (countSymbols)
```

Number of genes detected in sub-sample


```r
l2.sub1 <- data.frame(l2.sub1)
g2.sub1 <- rowsum(l2.sub1, genesymbols$symbol)
countSymbols.sub1 <- countSymbols(g2.sub1)
libs[colnames(l2.sub1),"genes.sub1"] <- (countSymbols.sub1)
```

### Table record
save the different tables produced for later analysis


```r
write.table(l2_NCki, "l2_NCki_1.txt", sep = "\t", quote=FALSE)
write.table(l2.sub1, "l2.sub1_NCki_1.txt", sep = "\t", quote=FALSE)
write.table(g2.sub1, 'g2.sub1_NCki_1.txt', sep="\t", quote=F)
write.table(libs, 'libs_NCki_1.txt', sep="\t", quote=F)
```
