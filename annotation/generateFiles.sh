#!/bin/bash

# Please look at https://gist.github.com/charles-plessy/9dbc8bc98fb773bf71b6
# For something way more up to date.

wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_14/gencode.v14.annotation.gtf.gz

zgrep -v '#' gencode.v14.annotation.gtf.gz |
  awk '{OFS="\t"} {$3 = "gene"} {print $1, $4, $5, $18, 0, $7}' |
  sed -e 's/[";]//g' > gencode.v14.annotation.genes.bed

GENCODE=gencode.v14.annotation.gtf.gz

zcat $GENCODE |
  grep -P '\tgene\t' |
  perl -F'\t' -anE 'say join "\t", @F[0], @F[3] - 1, @F[4], join ("_", @F[8] =~ /transcript_type "(.*?)"/, @F[8] =~ /transcript_name "(.*?)"/), 0, @F[6]' > annot.bed

zcat $GENCODE |
  grep -P '\texon\t' |
  perl -F'\t' -anE 'say join "\t", @F[0], @F[3] - 1, @F[4], "exon", 0, @F[6]' >> annot.bed

zcat $GENCODE |
  grep -P '\texon\t' |
  awk '{OFS="\t"}
    {print $1, $4 - 10, $4 + 10, "boundary", 0, $7}
    {print $1, $5 - 10, $5 + 10, "boundary", 0, $7}' >> annot.bed

zcat $GENCODE |
  grep -P '\tgene\t' |
  awk '{OFS="\t"} {
      if ($7=="+") {print $1,$4 - 100,$4 + 100,"TSS",0,$7}
      if ($7=="-") {print $1,$5 - 100,$5 + 100,"TSS",0,$7}
    }' >> annot.bed

zcat $GENCODE |
  grep -P '\tgene\t' |
  awk '{OFS="\t"} {
      if ($7=="+") {print $1,$5 - 100,$5 + 100,"TTS",0,$7}
      if ($7=="-") {print $1,$4 - 100,$4 + 100,"TTS",0,$7}
    }' >> annot.bed

md5sum -c md5sums.txt 
