origchrs = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', "X", "Y", 'GL000191.1', 'GL000192.1', 'GL000193.1', 'GL000194.1', 'GL000195.1', 'GL000196.1', 'GL000197.1', 'GL000198.1', 'GL000199.1', 'GL000200.1', 'GL000201.1', 'GL000202.1', 'GL000203.1', 'GL000204.1', 'GL000205.1', 'GL000206.1', 'GL000207.1', 'GL000208.1', 'GL000209.1', 'GL000210.1', 'GL000211.1', 'GL000212.1', 'GL000213.1', 'GL000214.1', 'GL000215.1', 'GL000216.1', 'GL000217.1', 'GL000218.1', 'GL000219.1', 'GL000220.1', 'GL000221.1', 'GL000222.1', 'GL000223.1', 'GL000224.1', 'GL000225.1', 'GL000226.1', 'GL000227.1', 'GL000228.1', 'GL000229.1', 'GL000230.1', 'GL000231.1', 'GL000232.1', 'GL000233.1', 'GL000234.1', 'GL000235.1', 'GL000236.1', 'GL000237.1', 'GL000238.1', 'GL000239.1', 'GL000240.1', 'GL000241.1', 'GL000242.1', 'GL000243.1', 'GL000244.1', 'GL000245.1', 'GL000246.1', 'GL000247.1', 'GL000248.1', 'GL000249.1']
newchrs1 = ["chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21", "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrX", "chrY"]
newchrs2 = ['chr1_gl000191_random', 'chr1_gl000192_random', 'chr4_gl000193_random', 'chr4_gl000194_random', 'chr7_gl000195_random', 'chr8_gl000196_random', 'chr8_gl000197_random', 'chr9_gl000198_random', 'chr9_gl000199_random', 'chr9_gl000200_random', 'chr9_gl000201_random', 'chr11_gl000202_random', 'chr17_gl000203_random', 'chr17_gl000204_random', 'chr17_gl000205_random', 'chr17_gl000206_random', 'chr18_gl000207_random', 'chr19_gl000208_random', 'chr19_gl000209_random', 'chr21_gl000210_random', 'chrUn_gl000211', 'chrUn_gl000212', 'chrUn_gl000213', 'chrUn_gl000214', 'chrUn_gl000215', 'chrUn_gl000216', 'chrUn_gl000217', 'chrUn_gl000218', 'chrUn_gl000219', 'chrUn_gl000220', 'chrUn_gl000221', 'chrUn_gl000222', 'chrUn_gl000223', 'chrUn_gl000224', 'chrUn_gl000225', 'chrUn_gl000226', 'chrUn_gl000227', 'chrUn_gl000228', 'chrUn_gl000229', 'chrUn_gl000230', 'chrUn_gl000231', 'chrUn_gl000232', 'chrUn_gl000233', 'chrUn_gl000234', 'chrUn_gl000235', 'chrUn_gl000236', 'chrUn_gl000237', 'chrUn_gl000238', 'chrUn_gl000239', 'chrUn_gl000240', 'chrUn_gl000241', 'chrUn_gl000242', 'chrUn_gl000243', 'chrUn_gl000244', 'chrUn_gl000245', 'chrUn_gl000246', 'chrUn_gl000247', 'chrUn_gl000248', 'chrUn_gl000249']
newchrs1.extend(newchrs2)
newchrs = newchrs1

newchrssave = newchrs
newchrs = origchrs
origchrs = newchrssave

inFI = open('/data/test1.sam')#'somatic-b37_Mutect2-WGS-panel-b37.vcf')
outFI = open('/data/test1_renamed.sam','w')#'somatic-hg19_Mutect2-WGS-panel-hg19.vcf','w')

linenum = 0
for line in inFI:
    linenum += 1
    if (linenum % 10000)==0:
        print(linenum)
    if linenum<100000:
        for j in range(len(origchrs)):
            if line.find("ID="+origchrs[j])!=-1:
                line = line.replace("ID="+origchrs[j],"ID="+newchrs[j])
    words = line.strip().split('\t')
    for j in range(len(origchrs)):
        line = line.replace(origchrs[j],newchrs[j]).strip()
        #if words[0]==origchrs[j]:
        #    words[0] = newchrs[j]
    #outFI.write('\t'.join(words)+'\n')
    if line.find("MT")==-1 and line.find("chrM")==-1:
        outFI.write(line+'\n')

inFI.close()
outFI.close()

#samtools view -S -b test_renamed.sam > test_renamed.bam
#samtools sort test_renamed.bam > test_renamed.sorted.bam
#samtools index test_renamed.sorted.bam

#wget https://storage.googleapis.com/genomics=public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta
#wget https://storage.googleapis.com/genomics=public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.fasta.fai
#~/hisat2-2.2.1/hisat2-build -p 16 Homo_sapiens_assembly38.fastq Homo_sapiens_assembly38
#mkdir Homo_sapiens_assembly38_hisat_index
#mv Homo_sapiens_assembly38*.ht2 Homo_sapiens_assembly38_hisat_index/
#cd Homo_sapiens_assembly38_hisat_index/
#~/hisat2-2.2.1/hisat2 -q -x Homo_sapiens_assembly38 -1 ../HM-baseline_TGACCA_L001_R1_001.fastq -2 ../HM-baseline_TGACCA_L001_R2_L001.fastq -S ../test_hg38.sam

#samtools view -S -b test_hg38.sam > test_hg38.bam
#samtools sort test_hg38.bam > test_hg38.sorted.bam
#samtools index test_hg38.sorted.bam

#~/gatk-4.2.6.1/gatk AddOrReplaceReadGroups I=test_hg38.sorted.bam O=test_hg38_rg.sorted.bam RGLB=lib1 RGPU=unit1 RGSM=20 RGID=4 RGPL=illumina
#samtools index test_hg38_rg.sorted.bam
#~/gatk-4.2.6.1/gatk PreprocessIntervals -R /data/Homo_sapiens_assembly38.fasta --padding 0 -imr OVERLAPPING_ONLY -O Homo_sapiens_assembly38.preprocessed.interval_list
#~/gatk-4.2.6.1/gatk CollectReadCounts -R /data/Homo_sapiens_assembly38.fasta -imr OVERLAPPING_ONLY --format TSV -L Homo_sapiens_assembly38.preprocessed.interval_list -I test_hg38_rg.sorted.bam -O test_hg38_rg.sorted.tsv
#~/gatk-4.2.6.1/gatk DetermineGermlineContigPloidy --model cohort-23wgs-20190213-contig-ploidy-model -I test_hg38_rg.sorted.tsv -O . --output-prefix ploidy-case --verbosity DEBUG
