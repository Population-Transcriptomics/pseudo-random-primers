Selection of Pseudo-random primers
======================================

### rRNA

Human:

hsu13369.fasta file produced with the command extractfeat -type rRNA U13369.gb. (from the EMBOSS package)

```{r}
HSU13369_3657_5527 [rRNA] Human ribosomal DNA complete repeating unit.
HSU13369_6623_6779 [rRNA] Human ribosomal DNA complete repeating unit.
HSU13369_7935_12969 [rRNA] Human ribosomal DNA complete repeating unit.
```

### Mitochondrial

Was extractfeat also used ? 

```{r}
NC_012920_648_1601 [rRNA] Homo sapiens mitochondrion, complete genome.
NC_012920_1671_3229 [rRNA] Homo sapiens mitochondrion, complete genome.
```

### Combination

```{r}
(cat nc_012920.fasta hsu13369.fasta | revseq -filter | grep -v '>' | perl -pe chomp ; echo) > ribo.txt 
```

R code
------

```{r}
acgt <- c('A', 'C', 'G', 'T')
LINKER <- 'CCCTATAAGATCGGAAGAGCGGTTCGGAGACCTTCAGTTCGACTA'
BARCODES <- scan('barcodes.txt', what='character')
RIBO <- scan('ribo.txt', what='character')  # See below in the wiki about the file 'ribo.txt'.
hexamers <- apply(expand.grid(acgt, acgt, acgt, acgt, acgt, acgt), 1, paste, collapse='')
hexamers <- data.frame(row.names=hexamers)
hexamers[,c('LINKER_0', 'LINKER_1', 'LINKER_2', 'LINKER_3', 'RIBO_0', 'RIBO_1', 'BARCODE')] <- c(rep(FALSE, 7))
hexamers[names(unlist(sapply(rownames(hexamers), function(X) {agrep(X, LINKER, 0, ignore.case=T)}))), "LINKER_0"] <- TRUE
hexamers[names(unlist(sapply(rownames(hexamers), function(X) {agrep(X, LINKER, 1, ignore.case=T)}))), "LINKER_1"] <- TRUE
hexamers[names(unlist(sapply(rownames(hexamers), function(X) {agrep(X, LINKER, 2, ignore.case=T)}))), "LINKER_2"] <- TRUE
hexamers[names(unlist(sapply(rownames(hexamers), function(X) {agrep(X, LINKER, 3, ignore.case=T)}))), "LINKER_3"] <- TRUE
hexamers[names(unlist(sapply(rownames(hexamers), function(X) {agrep(X, RIBO, 0, ignore.case=T)}))), "RIBO_0"] <- TRUE
hexamers[names(unlist(sapply(rownames(hexamers), function(X) {agrep(X, RIBO, 1, ignore.case=T)}))), "RIBO_1"] <- TRUE
hexamers[BARCODES, "BARCODE"] <- TRUE
```

```{r}
summary(hexamers)
LINKER_0        LINKER_1        LINKER_2       LINKER_3         RIBO_0         RIBO_1         BARCODE       
 Mode :logical   Mode :logical   Mode :logical   Mode:logical   Mode :logical   Mode:logical   Mode :logical  
 FALSE:4056      FALSE:3082      FALSE:259       TRUE:4096      FALSE:719       TRUE:4096      FALSE:4000     
 TRUE :40        TRUE :1014      TRUE :3837      NA's:0         TRUE :3377      NA's:0         TRUE :96       
 NA's :0         NA's :0         NA's :0                        NA's :0  
```

```{r}
with(hexamers, rownames(hexamers)[! (LINKER_2 | RIBO_0 | BARCODE)])
 [1] "GCCAAA" "AGCAAA" "AAACAA" "ACACAA" "TGCCAA" "CAAACA" "CACACA" "TGCACA"
 [9] "GTCACA" "TAGCCA" "GTGGCA" "TGTTTA" "ATTTTA" "CAAAAC" "CACAAC" "GCTAAC"
[17] "AACCAC" "CTACCC" "TACCCC" "CTAGCC" "CTGGCC" "TGTGCC" "ATTGCC" "CTACGC"
[25] "TATGGC" "TTGTGC" "ACCACG" "CACAGG" "ACTGTG" "TGCCAT" "TGGCAT" "GTGCAT"
[33] "TTGTAT" "ATTTAT" "TTTTAT" "TGGCGT" "TGTTGT" "ATTTGT" "TTGCTT" "TGTCTT"
```

Selection PS_Hb
======================================

### Haemoglobin sequences

alpha globin mRNA : http://www.ncbi.nlm.nih.gov/nuccore/NM_000558

beta globin mRNA : http://www.ncbi.nlm.nih.gov/nuccore/NM_000518

combine in one file Hb.txt


R Code
------

```{r}
acgt <- c('A', 'C', 'G', 'T')
Hb <- scan('Hb.txt', what='character')
hexamers <- apply(expand.grid(acgt, acgt, acgt, acgt, acgt, acgt), 1, paste, collapse='')
hexamers <- data.frame(row.names=hexamers)
hexamers[,c('Hb_0', 'Hb_1', 'Hb_2')] <- c(rep(FALSE,3 ))
hexamers[names(unlist(sapply(rownames(hexamers), function(X) {agrep(X, Hb, 0, ignore.case=T)}))), "Hb_0"] <- TRUE
hexamers[names(unlist(sapply(rownames(hexamers), function(X) {agrep(X, Hb, 1, ignore.case=T)}))), "Hb_1"] <- TRUE
hexamers[names(unlist(sapply(rownames(hexamers), function(X) {agrep(X, Hb, 2, ignore.case=T)}))), "Hb_2"] <- TRUE
```

```{r}
summary(hexamers)
    Hb_0            Hb_1           Hb_2        
 Mode :logical   Mode :logical   Mode:logical  
FALSE:3154      FALSE:33        TRUE:4096     
TRUE :942       TRUE :4063      NA's:0        
NA's :0         NA's :0                       
```

```{r}
with(hexamers, rownames(hexamers)[! (Hb_1)])
[1] "GTTAAA" "CGACAA" "GGATAA" "GTATAA" "CTACGA" "TATCGA" "CGAATA" "GATATA"
[9] "CGTATA" "GTACTA" "TACCTA" "ATCGTA" "CTCGTA" "TCGTTA" "TAAAAC" "TACAAC"
[17] "ATTTAC" "AAACCC" "TAATGC" "ATCTGC" "CTAATC" "ATTCCG" "CTATCG" "GATTCG"
[25] "TACGAT" "ATCGAT" "ATCTAT" "TCGTAT" "CTAATT" "TCCATT" "CCGATT" "TCGATT"
[33] "CGATTT"
```

Selection of 40N6 primers
======================================

R code
------

```{r}
acgt <- c('A', 'C', 'G', 'T')
hexamers <- apply(expand.grid(acgt, acgt, acgt, acgt, acgt, acgt), 1, paste, collapse='')
sample(hexamers,40)
 [1] "CCAGTC" "CCCTTC" "TTTTTT" "CTGTAC" "TGACCG" "TGTGAT" "AACCCT" "AGGCGG"
 [9] "TCGTCT" "CTACAA" "GTACGC" "CAGAAG" "GTGTCT" "GTGTGC" "AAGACT" "CGGGTA"
[17] "AAGAGA" "GAGGTG" "GCTCTT" "GGTGTG" "GCACGT" "TGAACT" "GGGGCG" "GAGAGG" 
[25] "CCTCAG" "TAAGTT" "ATCTGC" "ACTTAA" "CACAGC" "AGATGA" "GGTAGC" "AAGGCC" 
[33] "CGCAGG" "AACCTC" "CAGTTG" "ATTCCC" "AGATGG" "GCGGAC" "CTGGCG" "CTTCAC"
```