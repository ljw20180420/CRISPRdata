ext1up=${ext1up:-100}
ext1down=${ext1down:-27}
ext2up=${ext2up:-27}
ext2down=${ext2down:-100}
genome=${genome:-../genome/genome.fa}

getStrand()
{
    if [ $type = "del" ] || [ $type = "up" ] || [ $type = "sg1" ]
    then
        tstrand1=$strand1
    elif [ $type = "down" ] || [ $type = "dup" ]
    then
        tstrand1=$(sed "s/$strand1//" <<<"+-")
    else
        tstrand1=$strand2
    fi
    if [ $type = "del" ] || [ $type = "down" ] || [ $type = "sg2" ]
    then
        tstrand2=$strand2
    elif [ $type = "up" ] || [ $type = "dup" ]
    then
        tstrand2=$(sed "s/$strand2//" <<<"+-")
    else
        tstrand2=$strand1
    fi
}

getChrCut()
{
    if [ $type = "sg2" ]
    then
        tchr1=$chr2
        tcut1=$cut2
    else
        tchr1=$chr1
        tcut1=$cut1
    fi
    if [ $type = "sg1" ]
    then
        tchr2=$chr1
        tcut2=$cut1
    else
        tchr2=$chr2
        tcut2=$cut2
    fi
}

getRef()
{
    if [ $tstrand1 = "+" ]
    then
        start1=$(($tcut1 - $ext1up))
        end1=$(($tcut1 + $ext1down))
    else
        start1=$(($tcut1 - $ext1down))
        end1=$(($tcut1 + $ext1up))
    fi
    if [ $tstrand2 = "+" ]
    then
        start2=$(($tcut2 - $ext2up))
        end2=$(($tcut2 + $ext2down))
    else
        start2=$(($tcut2 - $ext2down))
        end2=$(($tcut2 + $ext2up))
    fi
    ref1=$(printf "%s\t%s\t%s\tref1\t.\t%s\n" $tchr1 $start1 $end1 $tstrand1 | bedtools getfasta -s -fi $genome -bed - | sed '1d')
    ref2=$(printf "%s\t%s\t%s\tref1\t.\t%s\n" $tchr2 $start2 $end2 $tstrand2 | bedtools getfasta -s -fi $genome -bed - | sed '1d')
    printf "0\t%s\t%s\t%s\t%s\t%s\n" $ref1 $ext1up $ext2up $ref2 ${#ref2}
}

CrickP5()
{
    temp=$tchr1
    tchr1=$tchr2
    tchr2=$temp
    temp=$tcut1
    tcut1=$tcut2
    tcut2=$temp
    temp=$tstrand1
    tstrand1=$tstrand2
    tstrand2=$temp
    tstrand1=$(sed "s/$tstrand1//" <<<"+-")
    tstrand2=$(sed "s/$tstrand2//" <<<"+-")
}

mkdir -p refs
while read name chr1 cut1 strand1 chr2 cut2 strand2
do
    for type in del up down dup sg1 sg2
    do
        getStrand
        getChrCut
        getRef >refs/$name-$type.Watson.ref
        CrickP5
        getRef >refs/$name-$type.Crick.ref
    done
done < ref_poses.tsv