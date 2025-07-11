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

getRefAndAlign()
{
    read fqFile refFile <<<$@
    read _ ref1 cut1 cut2 ref2 _ <${ref_path}/$refFile
    if [ ${ref1:$(($cut1 + 4)):2} == "GG" ]
    then
        PAM1="NGG"
    else
        PAM1="CCN"
    fi
    if [ ${ref2:$(($cut2 + 4)):2} == "GG" ]
    then
        PAM2="NGG"
    else
        PAM2="CCN"
    fi
    zcat ${data_path}/$fqFile.gz | sed -n '2~4p' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1, 0}' | rearrangement 3<${ref_path}/$refFile -s0 $s0 -s1 $s1 -s2 $s2 -u $u -v $v -ru $ru -rv $rv -qu $qu -qv $qv | gawk -f correct_micro_homology.awk -- ${ref_path}/$refFile $PAM1 $PAM2 | sed 'N;N;s/\n/\t/g' | sort -k2,2nr | gzip >${alg_path}/$fqFile.alg.gz
}
export -f getRefAndAlign

data_path=${DATA_DIR}/LE/data
ref_path=${DATA_DIR}/LE/refs
alg_path=${DATA_DIR}/LE/algs
mkdir -p ${alg_path}
parallel -a fq2ref.tsv --jobs 12 getRefAndAlign