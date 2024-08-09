#!/bin/bash

export s0=${s0:--6}
export s1=${s1:-4}
export s2=${s2:-2}
export u=${u:--3}
export v=${v:--9}
export ru=${ru:-0}
export rv=${rv:-0}
export qu=${qu:-0}
export qv=${qv:--5}
export minScoreTarget=${minScoreTarget:-30}
export minScorePair=${minScorePair:-100}
export minToMapShear=${minToMapShear:-30}

demultiplexAndAlign()
{
    read target refFile spliterTarget spliterPair <<<$@
    pair=${target//.R2/}
    rmDupFile=$(mktemp)
    removeDuplicates.sh data/$target data/$pair >$rmDupFile
    demultiplexFile=$(mktemp)
    spliterTarget=csvfiles/$spliterTarget spliterPair=csvfiles/$spliterPair minScoreTarget=$minScoreTarget minScorePair=$minScorePair demultiplex.sh $rmDupFile >$demultiplexFile
    sxCutR2AdapterFilterCumulate.sh $demultiplexFile $minToMapShear | rearrangement 3<refs/$refFile -s0 $s0 -s1 $s1 -s2 $s2 -u $u -v $v -ru $ru -rv $rv -qu $qu -qv $qv | gawk -f correct_micro_homology.awk -- refs/$refFile NGG NGG | sed 'N;N;s/\n/\t/g' | sort -k2,2nr | gzip >algs/${target%.gz}.alg.gz
}
export -f demultiplexAndAlign

parallel -a fq2ref.tsv --jobs 12 demultiplexAndAlign