#!/bin/bash

mkdir -p genome
# example: /home/abc/hg19/hg19.fa
read -ep "please input the full path of genome fasta file:" genome
ln -sf $genome ./genome/genome.fa
for ext in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2
do
    ln -sf ${genome%.*}.$ext ./genome/genome.$ext
done
chown $(id -un):$(id -gn) genome