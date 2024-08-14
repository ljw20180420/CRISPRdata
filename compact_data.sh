#!/bin/bash

concat_algs()
{
    for author in SX SJLJH LE
    do
        for alg in $(find data/$author/algs -type f)
        do
            bname=$(basename $alg)
            # ref1, ref2, cut1, cut2, author, file_name, ref1_end, ref2_start, random_insert, count
            zcat $alg | awk -F "\t" -v OFS="\t" -v author=$author -v file_name=${bname%.alg.gz} '{
                gsub("-", "", $18)
                ref1_len = match($18, /[acgtn]-*[acgtn]/)
                print toupper(substr($18, 1, ref1_len)), toupper(substr($18, ref1_len + 1)), $16, $17, author, file_name, $8, $11, $10, $2
            }'
        done
    done
}

concat_algs | sort -t $'\t' -k1,9 | awk -F "\t" '
function print_array1D(array1D, format)
{
    for (j in array1D)
    {
        printf(format, array1D[j])
        if (j < length(array1D))
            printf(",")
    }
}
function print_array2D(array2D, format)
{
    for (i in array2D)
    {
        printf("[")
        print_array1D(array2D[i], format)
        printf("]")
        if (i < length(array2D))
            printf(",")
    }
}
{
    new_ref = $1 != ref1 || $2 != ref2 || $3 != cut1 || $4 != cut2
    new_author = $5 != authors[au_len] || $6 != file_names[au_len]
    new_end = $7 != ref1_ends[au_len][end_len] || $8 != ref2_starts[au_len][end_len] || $9 != random_inserts[au_len][end_len]
    if (new_ref) {
        if (ref1) 
        {
            printf("{\"ref1\":\"%s\",\"ref2\":\"%s\",\"cut1\":%d,\"cut2\":%d,\"author\":[", ref1, ref2, cut1, cut2)
            print_array1D(authors, "\"%s\"")
            printf("],\"file_name\":[")
            print_array1D(file_names, "\"%s\"")
            printf("],\"ref1_end\":[")
            print_array2D(ref1_ends, "%d")
            printf("],\"ref2_start\":[")
            print_array2D(ref2_starts, "%d")
            printf("],\"random_insert\":[")
            print_array2D(random_inserts, "\"%s\"")
            printf("],\"count\":[")
            print_array2D(counts, "%d")
            printf("]}\n")
        }
        ref1 = $1
        ref2 = $2
        cut1 = $3
        cut2 = $4
        delete authors
        delete file_names
        delete ref1_ends
        delete ref2_starts
        delete random_inserts
        delete counts
        au_len = 0
    }
    if (new_ref || new_author)
    {
        ++au_len
        end_len = 0
        authors[au_len] = $5
        file_names[au_len] = $6
    }
    if (new_ref || new_author || new_end)
    {
        ++end_len
        ref1_ends[au_len][end_len] = $7
        ref2_starts[au_len][end_len] = $8
        random_inserts[au_len][end_len] = $9
    }
    counts[au_len][end_len] += $10
}
' | gzip >dataset.json.gz