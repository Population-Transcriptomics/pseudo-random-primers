

Analyze of the second experiment: NC12
======================================

Cluster and annotate in the shell (not in R)
--------------------------------------------


```bash
LIBRARY=NC12_1
BAMFILES=../Moirai/NC12_1.CAGEscan_short-reads.20150629125015/properly_paired_rmdup/*bam

level1.py -o /dev/stdout -f 66 -F 516 $BAMFILES |
  bgzip > $LIBRARY.l1.gz
cat <(zgrep \# -A1 $LIBRARY.l1.gz) <(zgrep -v \# $LIBRARY.l1.gz | sed '1d' |
  sort --field-separator $'\t' -k2.4,2n -k 2.4,2.4 -k3,3n -k4,4n -k5,5) |
  bgzip |
  sponge $LIBRARY.l1.gz
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
## Opening NC12_1.l1.gz
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

LIBRARY <- "NC12_1"
```

### Load data


```r
l2_NC12 <- read.osc(paste(LIBRARY,'l2','gz',sep='.'), drop.coord=T, drop.norm=T)
colnames(l2_NC12) <- sub('raw.NC12_1', 'NC12', colnames(l2_NC12))
colSums(l2_NC12)
```

```
##  NC12.HeLa_40N6_A  NC12.HeLa_40N6_B  NC12.HeLa_40N6_C    NC12.HeLa_PS_A    NC12.HeLa_PS_B    NC12.HeLa_PS_C 
##             12154             17411             20790             24065             27215             54835 
## NC12.HeLa_RanN6_A NC12.HeLa_RanN6_B NC12.HeLa_RanN6_C  NC12.THP1_40N6_A  NC12.THP1_40N6_B  NC12.THP1_40N6_C 
##             10944             35582             23215              9271             15299             15775 
##    NC12.THP1_PS_A    NC12.THP1_PS_B    NC12.THP1_PS_C NC12.THP1_RanN6_A NC12.THP1_RanN6_B NC12.THP1_RanN6_C 
##             21303             23454             37395             13356             58890             34922
```

### Normalization number of read per sample: l2.sub ; libs$genes.sub

In all the 3 libraries used, one contain only few reads tags. The smallest one has 8,708 counts. In order to make meaningful comparisons, all of them are subsapled to 8700 counts.


```r
set.seed(1)
l2.sub1 <- t(rrarefy(t(l2_NC12),min(8700)))
colSums(l2.sub1)
```

```
##  NC12.HeLa_40N6_A  NC12.HeLa_40N6_B  NC12.HeLa_40N6_C    NC12.HeLa_PS_A    NC12.HeLa_PS_B    NC12.HeLa_PS_C 
##              8700              8700              8700              8700              8700              8700 
## NC12.HeLa_RanN6_A NC12.HeLa_RanN6_B NC12.HeLa_RanN6_C  NC12.THP1_40N6_A  NC12.THP1_40N6_B  NC12.THP1_40N6_C 
##              8700              8700              8700              8700              8700              8700 
##    NC12.THP1_PS_A    NC12.THP1_PS_B    NC12.THP1_PS_C NC12.THP1_RanN6_A NC12.THP1_RanN6_B NC12.THP1_RanN6_C 
##              8700              8700              8700              8700              8700              8700
```

### Moirai statistics

Load the QC data produced by the Moirai workflow with which the libraries were processed. Sort in the same way as the l1 and l2 tables, to allow for easy addition of columns.


```r
libs <- loadMoiraiStats(multiplex = "NC12_1.multiplex.txt", summary = "../Moirai/NC12_1.CAGEscan_short-reads.20150629125015/text/summary.txt", pipeline = "CAGEscan_short-reads")
rownames(libs) <- sub('HeLa', 'NC12.HeLa', rownames(libs))
rownames(libs) <- sub('THP1', 'NC12.THP1', rownames(libs))
```

### Number of clusters

Count the number of unique L2 clusters per libraries after subsampling, and add
this to the QC table.  Each subsampling will give a different result, but the
mean result can be calculated by using the `rarefy` function at the same scale
as the subsampling.


```r
libs["l2.sub1"]     <- colSums(l2.sub1 > 0)
libs["l2.sub1.exp"] <- rarefy(t(l2_NC12), min(colSums(l2_NC12)))
```

### Richness

Richness should also be calculated on the whole data.


```r
libs["r100.l2"] <- rarefy(t(l2_NC12),100)
boxplot(data=libs, r100.l2 ~ group, ylim=c(92,100), las=1)
```

![](commandes_clean_files/figure-html/NC12_boxplot-1.png) 

### Hierarchical annotation

Differences of sampling will not bias distort the distribution of reads between annotations, so the non-subsampled library is used here.


```r
annot.l2 <- read.table(paste(LIBRARY,'l2','annot',sep='.'), head=F, col.names=c('id', 'feature'), row.names=1)
annot.l2 <- hierarchAnnot(annot.l2)

libs <- cbind(libs, t(rowsum(l2_NC12,  annot.l2[,'class']))) 
libs$samplename <- sub('HeLa', 'NC12_HeLa', libs$samplename)
libs$samplename <- sub('THP1', 'NC12_THP1', libs$samplename)
```

### Gene symbols used normalisation data


```r
genesymbols <- read.table(paste(LIBRARY,'l2','genes',sep='.'), col.names=c("cluster","symbol"), stringsAsFactors=FALSE)
rownames(genesymbols) <- genesymbols$cluster

g2 <- rowsum(l2_NC12, genesymbols$symbol)
countSymbols <- countSymbols(g2)

libs[colnames(l2_NC12),"genes"] <- (countSymbols)
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
write.table(l2_NC12, "l2_NC12_1.txt", sep = "\t", quote=FALSE)
write.table(l2.sub1, "l2.sub1_NC12_1.txt", sep = "\t", quote=FALSE)
write.table(g2.sub1, 'g2.sub1_NC12_1.txt', sep="\t", quote=F)
write.table(libs, 'libs_NC12_1.txt', sep="\t", quote=F)
```
