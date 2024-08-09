#!/bin/bash

ext1up=${ext1up:-100}
ext1down=${ext1down:-27}
ext2up=${ext2up:-27}
ext2down=${ext2down:-100}
genome=${genome:-$HOME/hg19_with_bowtie2_index/hg19.fa}

for sample in HPRT DCK
do
    for primer in U D
    do
        refFile=refs/$sample-$primer.ref
        >$refFile
        while read chr1 cut1 strand1 chr2 cut2 strand2 description
        do
            if ! grep -q $sample-$primer <<<$description
            then
                continue
            fi
            if [ $strand1 == "+" ]
            then
                start1=$(($cut1 - $ext1up))
                end1=$(($cut1 + $ext1down))
            else
                start1=$(($cut1 - $ext1down))
                end1=$(($cut1 + $ext1up))
            fi
            if [ $strand2 == "+" ]
            then
                start2=$(($cut2 - $ext2up))
                end2=$(($cut2 + $ext2down))
            else
                start2=$(($cut2 - $ext2down))
                end2=$(($cut2 + $ext2up))
            fi
            ref1=$(printf "%s\t%s\t%s\tref1\t.%s\n" $chr1 $start1 $end1 $strand1 | bedtools getfasta -fi $genome -bed - | sed '1d')
            ref2=$(printf "%s\t%s\t%s\tref1\t.%s\n" $chr2 $start2 $end2 $strand2 | bedtools getfasta -fi $genome -bed - | sed '1d')
            printf "0\t%s\t%s\t%s\t%s\t%s\n" $ref1 $ext1up $ext2up $ref2 ${#ref2} >>$refFile
        done <HPRTDCKUD.region
    done
done