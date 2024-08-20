#!/usr/bin/env -S gawk -f

function print_ref(      i)
{
    printf("{\"ref1\":\"%s\",\"ref2\":\"%s\",\"cuts\":[", ref["ref1"], ref["ref2"])
    for (i=1; i<=length(ref["cuts"]); ++i)
    {
        print_cut(ref["cuts"][i])
        if (i < length(ref["cuts"]))
            printf(",")
    }
    printf("]}\n")
}

function print_cut(cut,      i)
{
    printf("{\"cut1\":%d,\"cut2\":%d,\"authors\":[", cut["cut1"], cut["cut2"])
    for (i=1; i<=length(cut["authors"]); ++i)
    {
        print_author(cut["authors"][i])
        if (i < length(cut["authors"]))
            printf(",")
    }
    printf("]}")
}

function print_author(author,      i)
{
    printf("{\"author\":\"%s\",\"files\":[", author["author"])
    for (i=1; i<=length(author["files"]); ++i)
    {
        print_file(author["files"][i])
        if (i < length(author["files"]))
            printf(",")
    }
    printf("]}")
}

function print_file(file,      i)
{
    printf("{\"file\":\"%s\",", file["file"])
    printf("\"ref1_end\":[")
    for (i=1; i<=length(file["ref1_ends"]); ++i)
    {
        printf("%d", file["ref1_ends"][i])
        if (i < length(file["ref1_ends"]))
            printf(",")
    }
    printf("],\"ref2_start\":[")
    for (i=1; i<=length(file["ref2_starts"]); ++i)
    {
        printf("%d", file["ref2_starts"][i])
        if (i < length(file["ref2_starts"]))
            printf(",")
    }
    printf("],\"random_insert\":[")
    for (i=1; i<=length(file["random_inserts"]); ++i)
    {
        printf("\"%s\"", file["random_inserts"][i])
        if (i < length(file["random_inserts"]))
            printf(",")
    }
    printf("],\"count\":[")
    for (i=1; i<=length(file["counts"]); ++i)
    {
        printf("%d", file["counts"][i])
        if (i < length(file["counts"]))
            printf(",")
    }
    printf("]}")
}

function process_ref()
{
    new_ref = $1 != ref["ref1"] || $2 != ref["ref2"]
    if (new_ref)
    {
        if (ref["ref1"])
            print_ref()
        delete ref
        ref["ref1"] = $1
        ref["ref2"] = $2
        # declare and delete to initialize an empty array
        ref["cuts"][1]
        delete ref["cuts"][1]
    }
    process_cut(ref["cuts"])
}

function process_cut(cuts,      i)
{
    i = length(cuts)
    new_cut = new_ref || $3 != cuts[i]["cut1"] || $4 != cuts[i]["cut2"]
    if (new_cut)
    {
        ++i
        cuts[i]["cut1"] = $3
        cuts[i]["cut2"] = $4
        # declare and delete to initialize an empty array
        cuts[i]["authors"][1]
        delete cuts[i]["authors"][1]
    }
    process_author(cuts[i]["authors"])
}

function process_author(authors,      i)
{
    i = length(authors)
    new_author = new_cut || $5 != authors[i]["author"] 
    if (new_author)
    {
        ++i
        authors[i]["author"] = $5
        # declare and delete to initialize an empty array
        authors[i]["files"][1]
        delete authors[i]["files"][1]
    }
    process_file(authors[i]["files"])
}

function process_file(files,      i)
{
    i = length(files)
    new_file = new_author || $6 != files[i]["file"]
    if (new_file)
    {
        ++i
        files[i]["file"] = $6
        # declare and delete to initialize an empty array
        files[i]["ref1_ends"][1]
        delete files[i]["ref1_ends"][1]
        files[i]["ref2_starts"][1]
        delete files[i]["ref2_starts"][1]
        files[i]["random_inserts"][1]
        delete files[i]["random_inserts"][1]
        files[i]["counts"][1]
        delete files[i]["counts"][1]
    }
    process_edit(files[i]["ref1_ends"], files[i]["ref2_starts"], files[i]["random_inserts"], files[i]["counts"])
}

function process_edit(ref1_ends, ref2_starts, random_inserts, counts,      i)
{
    i = length(ref1_ends)
    new_edit = new_file || $7 != ref1_ends[i] || $8 != ref2_starts[i] || $9 != random_inserts[i]
    if (new_edit)
    {
        ++i
        ref1_ends[i] = $7
        ref2_starts[i] = $8
        random_inserts[i] = $9
    }
    counts[i] += $10
}

BEGIN{
    FS="\t"
}
{
    if ($11 < min_score)
        next
    process_ref()
}
END{
    print_ref()
}