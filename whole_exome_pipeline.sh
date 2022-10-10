wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Homo_sapiens_assembly19.fasta Homo_sapiens_assembly19.fasta
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Homo_sapiens_assembly19.fasta.fai Homo_sapiens_assembly19.fasta.fai
~/hisat2-2.2.1/hisat2-build -p 16 Homo_sapiens_assembly19.fastq Homo_sapiens_assembly19
mkdir Homo_sapiens_assembly19_hisat_index
mv Homo_sapiens_assembly19*.ht2 Homo_sapiens_assembly19_hisat_index/
cd Homo_sapiens_assembly19_hisat_index/
~/hisat2-2.2.1/hisat2 -q -x Homo_sapiens_assembly19 -1 ../HM-baseline_TGACCA_R1.fastq -2 ../HM-baseline_TGACCA_R2.fastq | samtools view -S -b > HM-baseline_TGACCA.bam

samtools sort HM-baseline_TGACCA.bam > HM-baseline_TGACCA.sorted.bam
~/gatk-4.2.6.1/gatk AddOrReplaceReadGroups I=HM-baseline_TGACCA.sorted.bam O=HM-baseline_TGACCA_rg.sorted.bam RGLB=lib1 RGPU=unit1 RGSM=20 RGID=4 RGPL=illumina
samtools index HM-baseline_TGACCA_rg.sorted.bam

python3 -i splitChromosomes.py

bigBedToBed < /data/KAPA_HyperExome_hg19_capture_targets.bb > /data/KAPA_HyperExome_hg19_capture_targets.bed

~/FREEC-11.6/src/freec -conf ~/FREEC-11.6/data/config_exome_HM-baseline.txt
