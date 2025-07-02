#!/bin/bash

# 切换运行路径到脚本路径
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

mkdir -p rnafold
for file in $(find csvfiles -name "*.csv")
do
    cut -d, -f2 \
        < $file |
    cut -c 21-40 |
    sed -r 's/(.+)/>\1\n\1/' |
    RNAfold --noPS \
        > rnafold/$(basename $file).rnafold.sgRNA
    cut -d, -f2 \
        < $file |
    cut -c 21- |
    sed -r 's/[ACGTN]+\s+$//' |
    sed -r 's/^([ACGTN]{20})/>\1\n/' |
    RNAfold --noPS \
        > rnafold/$(basename $file).rnafold.sgRNA+scaffold
done