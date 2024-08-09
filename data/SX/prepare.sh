#!/bin/bash

ext1up=${ext1up:-100}
ext1down=${ext1down:-27}
ext2up=${ext2up:-27}
ext2down=${ext2down:-100}
genome=${genome:-../genome/genome.fa}
bowtie2index=${bowtie2index:-$HOME/hg19_with_bowtie2_index/hg19}

mkdir -p refs
for csvfile in $(ls csvfiles/*.csv)
do
    getSxCsvFileRef.sh $csvfile $genome $bowtie2index $ext1up $ext1down $ext2up $ext2down >refs/$(basename $csvfile).ref
    sxExtractSpliter.sh $csvfile >$csvfile.target.fa 3>$csvfile.pair.fa
    bowtie2-build $csvfile.target.fa $csvfile.target.fa
    bowtie2-build $csvfile.pair.fa $csvfile.pair.fa
done

>fq2ref.tsv
for target in $(ls data/*.R2.fq.gz)
do
    GAN=$(cut -d- -f2 <<<$target | head -c2 | dd conv=ucase 2>/dev/null)
    case $GAN in
        G?)
            spliterTarget="final_hgsgrna_libb_all_0811_NGG_scaffold_nor_${GAN}.csv.target.fa"
            spliterPair="final_hgsgrna_libb_all_0811_NGG_scaffold_nor_${GAN}.csv.pair.fa"
            refFile="final_hgsgrna_libb_all_0811_NGG_scaffold_nor_${GAN}.csv.ref"
            ;;
        A?)
            spliterTarget="final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_${GAN}.csv.target.fa"
            spliterPair="final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_${GAN}.csv.pair.fa"
            refFile="final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_${GAN}.csv.ref"
            ;;
    esac
    printf "%s\t%s\t%s\t%s\n" $(basename $target) $refFile $spliterTarget $spliterPair >>fq2ref.tsv
done