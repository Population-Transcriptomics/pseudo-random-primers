# Selection of Pseudo-random primers



Selection of Pseudo-random primers against rRNA
-----------------------------------------------

### Human rRNA reference sequences

We extracted the human rRNA sequences from the reference locus
[U13369](http://www.ncbi.nlm.nih.gov/nuccore/U13369) in a file called
`hsu13369.fasta`, produced with the command `extractfeat -type rRNA U13369.gb`
(from the [EMBOSS](http://emboss.sourceforge.net/) package).

```
HSU13369_3657_5527 [rRNA] Human ribosomal DNA complete repeating unit.
HSU13369_6623_6779 [rRNA] Human ribosomal DNA complete repeating unit.
HSU13369_7935_12969 [rRNA] Human ribosomal DNA complete repeating unit.
```

### Mitochondrial rRNA

We extracted human mitochondrial rRNA sequences in the same way, from the
reference locus [NC_012920](http://www.ncbi.nlm.nih.gov/nuccore/NC_012920).

```
NC_012920_648_1601 [rRNA] Homo sapiens mitochondrion, complete genome.
NC_012920_1671_3229 [rRNA] Homo sapiens mitochondrion, complete genome.
```

### Combination

The reference sequences were then reverse-complemented and combined in a
plain text file with the following command.

```
(cat nc_012920.fasta hsu13369.fasta | revseq -filter | grep -v '>' | perl -pe chomp ; echo) > ribo.txt 
```

### Selection in `R`

Reverse-complement of the linker sequence from
[Harbers _et al._, 2013](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-665),
but without barcode and fingerprints.


```r
LINKER <- 'CCCTATAAGATCGGAAGAGCGGTTCGGAGACCTTCAGTTCGACTA'
```

Barcode sequences from Poulain _et al._, 2016 (in press).  Note that when barcode
sequences are introduced after the reverse-transcription, it is not necessary to filter
them out.


```r
BARCODES <- scan('barcodes.txt', what='character')
```

rRNA sequence (see above for details)


```r
RIBO <- scan('ribo.txt', what='character')
```

Creation of a table representing every possible hexamer.


```r
acgt <- c('A', 'C', 'G', 'T')
hexamers <- apply(expand.grid(acgt, acgt, acgt, acgt, acgt, acgt), 1, paste, collapse='')
hexamers <- data.frame(row.names=hexamers)
```


```r
hexamers[,c('LINKER_0', 'LINKER_1', 'LINKER_2', 'LINKER_3', 'RIBO_0', 'RIBO_1', 'BARCODE')] <- c(rep(FALSE, 7))

matchHexamers <- function(reference, mismatches) {
  names( unlist( sapply( rownames(hexamers)
                       , function(X) {agrep( X
                                           , reference
                                           , mismatches
                                           , ignore.case = TRUE)}
                       )
               )
        )
}

hexamers[matchHexamers(LINKER, 0), "LINKER_0"] <- TRUE
hexamers[matchHexamers(LINKER, 1), "LINKER_1"] <- TRUE
hexamers[matchHexamers(LINKER, 2), "LINKER_2"] <- TRUE
hexamers[matchHexamers(LINKER, 3), "LINKER_3"] <- TRUE
hexamers[matchHexamers(RIBO, 0),   "RIBO_0"]   <- TRUE
hexamers[matchHexamers(RIBO, 1),   "RIBO_1"]   <- TRUE
hexamers[BARCODES,                 "BARCODE"]  <- TRUE
```

Here, barcodes are filtered out by matching row names directly, since in our nanoCAGE
method they are hexamers.  Note that the `matchHexamers` function does not expect line
breaks, so passing the a list of longer barcodes (like Illumina/Nextera indexes) will
not produce the expected results.


```r
summary(hexamers)
```

```
##   LINKER_0        LINKER_1        LINKER_2       LINKER_3         RIBO_0         RIBO_1       
##  Mode :logical   Mode :logical   Mode :logical   Mode:logical   Mode :logical   Mode:logical  
##  FALSE:4056      FALSE:3082      FALSE:259       TRUE:4096      FALSE:719       TRUE:4096     
##  TRUE :40        TRUE :1014      TRUE :3837      NA's:0         TRUE :3377      NA's:0        
##  NA's :0         NA's :0         NA's :0                        NA's :0                       
##   BARCODE       
##  Mode :logical  
##  FALSE:4000     
##  TRUE :96       
##  NA's :0
```


```r
with(hexamers, rownames(hexamers)[! (LINKER_2 | RIBO_0 | BARCODE)])
```

```
##  [1] "GCCAAA" "AGCAAA" "AAACAA" "ACACAA" "TGCCAA" "CAAACA" "CACACA" "TGCACA" "GTCACA" "TAGCCA"
## [11] "GTGGCA" "TGTTTA" "ATTTTA" "CAAAAC" "CACAAC" "GCTAAC" "AACCAC" "CTACCC" "TACCCC" "CTAGCC"
## [21] "CTGGCC" "TGTGCC" "ATTGCC" "CTACGC" "TATGGC" "TTGTGC" "ACCACG" "CACAGG" "ACTGTG" "TGCCAT"
## [31] "TGGCAT" "GTGCAT" "TTGTAT" "ATTTAT" "TTTTAT" "TGGCGT" "TGTTGT" "ATTTGT" "TTGCTT" "TGTCTT"
```


Selection of PS primers against hemoglobin
------------------------------------------

### Hemoglobin sequences

 - alpha globin mRNA: <https://www.ncbi.nlm.nih.gov/nuccore/NM_000558>.
 - beta globin mRNA: <https://www.ncbi.nlm.nih.gov/nuccore/NM_000518>.

Combined in one file `Hb.txt` (without FASTA headers).

### R Code


```r
acgt <- c('A', 'C', 'G', 'T')
Hb <- scan('Hb.txt', what='character')
hexamers <- apply(expand.grid(acgt, acgt, acgt, acgt, acgt, acgt), 1, paste, collapse='')
hexamers <- data.frame(row.names=hexamers)
hexamers[,c('Hb_0', 'Hb_1', 'Hb_2')] <- c(rep(FALSE,3 ))

hexamers[matchHexamers(Hb, 0), "Hb_0"] <- TRUE
hexamers[matchHexamers(Hb, 1), "Hb_1"] <- TRUE
hexamers[matchHexamers(Hb, 2), "Hb_2"] <- TRUE
```


```r
summary(hexamers)
```

```
##     Hb_0            Hb_1           Hb_2        
##  Mode :logical   Mode :logical   Mode:logical  
##  FALSE:3154      FALSE:33        TRUE:4096     
##  TRUE :942       TRUE :4063      NA's:0        
##  NA's :0         NA's :0
```


```r
with(hexamers, rownames(hexamers)[! (Hb_1)])
```

```
##  [1] "GTTAAA" "CGACAA" "GGATAA" "GTATAA" "CTACGA" "TATCGA" "CGAATA" "GATATA" "CGTATA" "GTACTA"
## [11] "TACCTA" "ATCGTA" "CTCGTA" "TCGTTA" "TAAAAC" "TACAAC" "ATTTAC" "AAACCC" "TAATGC" "ATCTGC"
## [21] "CTAATC" "ATTCCG" "CTATCG" "GATTCG" "TACGAT" "ATCGAT" "ATCTAT" "TCGTAT" "CTAATT" "TCCATT"
## [31] "CCGATT" "TCGATT" "CGATTT"
```


Selection of 40 random hexamers
-------------------------------


```r
acgt <- c('A', 'C', 'G', 'T')
hexamers <- apply(expand.grid(acgt, acgt, acgt, acgt, acgt, acgt), 1, paste, collapse='')
set.seed(1)
sample(hexamers,40)
```

```
##  [1] "TTTAAC" "TATTCC" "CGGACG" "CCAGGT" "CGTATA" "TGCCGT" "TCCATT" "GTAGGG" "TGAAGG" "ATTTAA"
## [11] "CGACTA" "CACTGA" "CCTTGG" "AAGAGC" "GCACAT" "TTGTTC" "TTGCTG" "CTATTT" "CTAAGC" "CAGCAT"
## [21] "CAGTGT" "AAGCTA" "GTCCGG" "TTTTCA" "AAACAC" "TAGAGC" "GCTAAA" "TACAGC" "CACTCT" "AGGCCC"
## [31] "AGGGTC" "CCAGCG" "CCCTTC" "ACTTGA" "AAGACT" "GGCGGG" "AGCGAT" "GCTGCA" "AGTCTG" "ACAGGC"
```

```r
rm(.Random.seed)
```

In our paper, the result of the random selection was:

```
 [1] "CCAGTC" "CCCTTC" "TTTTTT" "CTGTAC" "TGACCG" "TGTGAT" "AACCCT" "AGGCGG"
 [9] "TCGTCT" "CTACAA" "GTACGC" "CAGAAG" "GTGTCT" "GTGTGC" "AAGACT" "CGGGTA"
[17] "AAGAGA" "GAGGTG" "GCTCTT" "GGTGTG" "GCACGT" "TGAACT" "GGGGCG" "GAGAGG" 
[25] "CCTCAG" "TAAGTT" "ATCTGC" "ACTTAA" "CACAGC" "AGATGA" "GGTAGC" "AAGGCC" 
[33] "CGCAGG" "AACCTC" "CAGTTG" "ATTCCC" "AGATGG" "GCGGAC" "CTGGCG" "CTTCAC"
```
