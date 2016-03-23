#!/bin/sh

getAndUnpack() {
  wget https://zenodo.org/record/48114/files/$1
  tar xvfz $1
  rm $1
}

getAndUnpack NC12_1.CAGEscan_short-reads.20150629125015.tar.gz
getAndUnpack NC16-17_1.CAGEscan_short-reads.20150625154740.tar.gz
getAndUnpack NC22b.CAGEscan_short-reads.20150625152335.tar.gz
getAndUnpack NCms10058_1.CAGEscan_short-reads.20150625154711.tar.gz
