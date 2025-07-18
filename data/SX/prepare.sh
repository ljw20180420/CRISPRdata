#!/bin/bash

ext1up=${ext1up:-100}
ext1down=${ext1down:-27}
ext2up=${ext2up:-27}
ext2down=${ext2down:-100}

data_path=${DATA_DIR}/SX/data 
ref_path=${DATA_DIR}/SX/refs
mkdir -p ${ref_path}
for csvfile in $(ls csvfiles/*.csv)
do
    getSxCsvFileRef.md $csvfile ${GENOME} ${BOWTIE2_INDEX} $ext1up $ext1down $ext2up $ext2down >${ref_path}/$(basename $csvfile).ref
    sxExtractSpliter.md $csvfile >$csvfile.target.fa 3>$csvfile.pair.fa
    bowtie2-build $csvfile.target.fa $csvfile.target.fa
    bowtie2-build $csvfile.pair.fa $csvfile.pair.fa
done

>fq2ref.tsv
for target in $(ls ${data_path}/*.R2.fq.gz)
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