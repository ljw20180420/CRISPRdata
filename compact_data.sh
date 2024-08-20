#!/bin/bash

concat_algs()
{
    for author in SX SJLJH LE
    do
        for alg in $(find data/$author/algs -type f)
        do
            bname=$(basename $alg)
            # ref1, ref2, cut1, cut2, author, file, ref1_end, ref2_start, random_insert, count, score
            zcat $alg | awk -F "\t" -v OFS="\t" -v author=$author -v file=${bname%.alg.gz} '{
                gsub("-", "", $18)
                ref1_len = match($18, /[acgtn]-*[acgtn]/)
                print toupper(substr($18, 1, ref1_len)), toupper(substr($18, ref1_len + 1)), $16, $17 - ref1_len, author, file, $8, $11 - ref1_len, toupper($10), $2, $3
            }'
        done
    done
}

concat_algs | sort --parallel 24 -t $'\t' -k1,9 | ./to_json.awk -v min_score=0 | gzip >dataset.json.gz