import bioframe
import requests

### download 
# file='https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/rmsk.txt.gz'
# r = requests.get(file)
# with open(os.path.basename(file), 'wb') as f:
#     f.write(r.content)

### read map analysis
genome_file = '/home/ljw/hg19_with_bowtie2_index/hg19.fa'
excu = '/home/ljw/wuqiang/lierlib/Rearrangement/build/rearrangement'
datadir = '/home/ljw/wuqiang/sxdata'
fastq_files_dict = {}
fastq_files = []
for sample in ['DCK','HPRT']:
    for primer in ['U','D']:
        for Ri in ['R1', 'R2']:
            fastq_files_dict[f'{sample}_{primer}_{Ri}'] = [os.path.join(datadir,file) for file in os.listdir(datadir) if os.path.isfile(os.path.join(datadir,file)) and file.find(f'-{sample}-{primer}-')!=-1 and file.endswith(f".{Ri}.fastq.gz")]
            fastq_files += fastq_files_dict[f'{sample}_{primer}_{Ri}']

rmsk = bioframe.read_table('rmsk.txt.gz')[[5, 6, 7]]
rmsk.columns = ['chrom', 'start', 'end']

cts = [{'sample' : 'HPRT', 'chr' : 'chrX', 'sgrnaU' : "GTAAGTAAGATCTTAAAATGAGG", 'sgrnaD' : "CCTGATTTTATTTCTGTAGGACT"}, {'sample' : 'DCK', 'chr' : 'chr4', 'sgrnaU' : "ATATTTAGAACTCTTTTCAGTGG", 'sgrnaD' : "CCCTGCCTTTTTCTTCCATCTCT"}]

genome = bioframe.load_fasta(genome_file)
for ct in cts:
    chrgenome = genome[ct['chr']].ff.fetch(ct['chr'])
    ct['cutU'] = chrgenome.find(ct['sgrnaU'])+len(ct['sgrnaU'])-6
    ct['cutD'] = chrgenome.find(ct['sgrnaD'])+6
    ct['rmsk'] = bioframe.select(rmsk, (ct['chr'], ct['cutU']-3000, ct['cutD']+3000))
    # ct['crmsk'] = pandas.DataFrame({'chrom' : [ct['chr']] * (len(ct['rmsk'])+1), 'start' : [ct['cutU']-3000] + ct['rmsk']['end'].values.tolist(), 'end' : ct['rmsk']['start'].values.tolist() + [ct['cutD']+3000]})

threads = 24
inner, outter = 100000, 110000
minseg = 20
absTini = 10
ext_file = [os.path.join(datadir,file) for file in os.listdir(datadir) if os.path.isfile(os.path.join(datadir,file)) and file.endswith(f".fastq.gz.fa.ext")]
for fastq_file in fastq_files:
    if fastq_file.endswith(".R2.fastq.gz"):
        continue
    
    # if f"{fastq_file}.fa.ext" in ext_file:
    #     continue

    for ct in cts:
        if fastq_file.find(ct['sample'])!=-1:
            sample = ct['sample']
            chrs1, chrs2, cutU, cutD = [ct['chr']]*3, [ct['chr']]*3, ct['cutU'], ct['cutD']
            refext = ct['cutD'] - ct['cutU'] + 150
            ctt = ct
            break
    
    if fastq_file.find(f'-{sample}-U-')!=-1:
        cuts1, cuts2, strands1, strands2 = [cutU]*3, [cutU,cutD,cutD], ['+','+','+'], ['+','+','-']
    else:
        cuts1, cuts2, strands1, strands2 = [cutD]*3, [cutD,cutU,cutU], ['-','-','-'], ['-','-','+']

    # simple check
    if not (chrs1[0]==chrs1[1] and chrs1[1]==chrs1[2]):
        raise Exception("chrs1 must all be the same")
    cut12 = list(set(cuts1 + cuts2))
    if len(cut12)!=2:
        raise Exception("total cut positions must be two")

    # unify reads
    # count_file, fasta_file = count_fastq(fastq_file)
    count_file, fasta_file = f'{fastq_file}.txt', f'{fastq_file}.fa'
    best_file = f'{fasta_file}.indel'