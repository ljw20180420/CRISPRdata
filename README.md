# Introduction
This is a repository use rearr https://github.com/ljw20180420/sx_lcy to align CRISPR data chimerically.

# Select genome
If multiple sgRNAs targeting various genome sites are used, demultiplex will be necessary. Demultiplex needs the genome bowtie2 index. Select the genome by running the following command in the project directory.
```bash
./select_genome.sh
```

# Alignment
We align data in `data/LE/` as follows.
```bash
cd data/LE
./prepare.sh
./align.sh
cd -
```
Those in `data/SJLJH` and `data/SX` are aligned similarly. `./prepare.sh` is mainly used to generate references around CRISPR targeting sites.

# Collect alignments as a JSON file
We collect alignments in `data/LE/algs`, `data/SJLJH/algs` and `data/SX/algs` by
```bash
./compact_data.sh
```
This generates a small `dataset.json.gz` file. The huggingface dataset for CRISPR AI training only depends on `dataset.json.gz`.
