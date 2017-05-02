# 2.1) Correct HCT116 somatic mutations on HG38

cat P007_48_DMSO_DHG05826_HNHCCCCXX_L3_1.fq.gz P007_48_DMSO_DHG05826_HNHCCCCXX_L4_1.fq.gz P007_48_DMSO_DHG05826_HNHCCCCXX_L5_1.fq.gz P007_48_DMSO_DHG05826_HNHCCCCXX_L6_1.fq.gz P007_48_DMSO_DHG05826_HNHCCCCXX_L7_1.fq.gz P007_48_DMSO_DHG05826_HNHCCCCXX_L8_1.fq.gz > P007_48_DMSO_1.fq.gz

cat P007_48_DMSO_DHG05826_HNHCCCCXX_L3_2.fq.gz P007_48_DMSO_DHG05826_HNHCCCCXX_L4_2.fq.gz P007_48_DMSO_DHG05826_HNHCCCCXX_L5_2.fq.gz P007_48_DMSO_DHG05826_HNHCCCCXX_L6_2.fq.gz P007_48_DMSO_DHG05826_HNHCCCCXX_L7_2.fq.gz P007_48_DMSO_DHG05826_HNHCCCCXX_L8_2.fq.gz > P007_48_DMSO_2.fq.gz

###
##

#

mkdir 2.1_mutation_correction
cd 2.1_mutation_correction

# 2.1.1) Variant Calling

# BWA installation
git clone https://github.com/lh3/bwa.git
cd bwa
ln -s /home/rtm/myPrograms/bwa-0.7.15/bwa /home/rtm/bin/bwa

# BWA generating genome index
bwa index -a bwtsw ~/resources/hg38/star/genome.fa

# BWA mapping
bwa mem -M -t 16 ~/resources/hg38/star/genome.fa /home/rtm/CY/RawData/P007_48_DMSO/P007_48_DMSO_1.fq.gz /home/rtm/CY/RawData/P007_48_DMSO/P007_48_DMSO_2.fq.gz > HCT116_DMSO_48h.sam
# -M ensures compatibility with picard
# -t number of threads

# SAM TO BAM + SORT
samtools view -Sb HCT116_DMSO_48h.sam |samtools sort - HCT116_DMSO_48h

# In order to reduce the number of miscalls of INDELs in your data it is helpful to realign your raw gapped alignment with the Broad’s GATK Realigner.

# broadBundle
 wget -r 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/'
# Fasta index
samtools faidx genome.fa
# Creating genome dir
java -jar ~/myPrograms/picard/build/libs/picard.jar CreateSequenceDictionary R= genome.fa O= genome.dict
# Add RG
java -jar ~/myPrograms/picard/build/libs/picard.jar AddOrReplaceReadGroups INPUT=HCT116_DMSO_48h.bam OUTPUT=HCT116_DMSO_48h_addRG.bam RGID=HNHCCCCXX RGLB= Merged RGPL=illumina RGPU=HNHCCCCXX RGSM=sample1 
# Index
samtools index HCT116_DMSO_48h_addRG.bam

# ReMapping
java -Xmx50g -jar ~/myPrograms/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ~/resources/hg38/star/genome.fa \
-I HCT116_DMSO_48h_addRG.bam -o lane.intervals \
--known ~/resources/ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf

# INDEL ReMapping
java -Xmx60g -jar  ~/myPrograms/GenomeAnalysisTK.jar -T IndelRealigner -nct 14 -R ~/resources/hg38/star/genome.fa \
-I HCT116_DMSO_48h_addRG.bam -targetIntervals lane.intervals \
-known ~/resources/ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf -o HCT116_DMSO_48h_addRG_realigned.bam

#Base recallibration 
java -Xmx60g -jar ~/myPrograms/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 14 -R ~/resources/hg38/star/genome.fa -knownSites ~/resources/ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf -I HCT116_DMSO_48h_addRG_realigned.bam -o HCT116_DMSO_48h_recal.table
java -Xmx60g -jar ~/myPrograms/GenomeAnalysisTK.jar -T PrintReads -nct 14 -R ~/resources/hg38/star/genome.fa -I HCT116_DMSO_48h_addRG_realigned.bam --BQSR HCT116_DMSO_48h_recal.table -o  HCT116_DMSO_48h_addRG_realigned_recalibrated.bam

# Remove duplicates
java -Xmx60g -jar ~/myPrograms/picard/build/libs/picard.jar MarkDuplicates  VALIDATION_STRINGENCY=LENIENT M=MarkDup_metrics.txt INPUT=HCT116_DMSO_48h_addRG_realigned_recalibrated.bam OUTPUT=HCT116_DMSO_48h_addRG_realigned_recalibrated_rmdup.bam

#INDEX
samtools index HCT116_DMSO_48h_addRG_realigned_recalibrated_rmdup.bam

# Variant calling pileup
~/myPrograms/samtools-1.3.1/samtools mpileup -ugf ~/resources/hg38/star/genome.fa HCT116_DMSO_48h_addRG_realigned_recalibrated_rmdup.bam | bcftools call -vmO z -o HCT116_DMSO_48h.vcf.gz

# VCF index
tabix -p vcf HCT116_DMSO_48h.vcf.gz

# Extract SNP and INDELS
java -jar ~/myPrograms/GenomeAnalysisTK.jar -T SelectVariants \
-R ~/resources/hg38/star/genome.fa \
-V HCT116_DMSO_48h.vcf.gz \
-selectType SNP \
-o HCT116_DMSO_48h_SNP.vcf.gz 

java -jar ~/myPrograms/GenomeAnalysisTK.jar -T SelectVariants \
-R ~/resources/hg38/star/genome.fa \
-V HCT116_DMSO_48h.vcf.gz \
-selectType INDEL \
-o HCT116_DMSO_48h_INDEL.vcf.gz 

# Filtering

java -jar ~/myPrograms/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ~/resources/hg38/star/genome.fa \
-V HCT116_DMSO_48h_SNP.vcf.gz \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filterName "snp_filter" \
-o HCT116_DMSO_48h_SNP_filtered.vcf.gz 

java -jar ~/myPrograms/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ~/resources/hg38/star/genome.fa \
-V HCT116_DMSO_48h_INDEL.vcf.gz \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filterName "indel_filter" \
-o HCT116_DMSO_48h_INDEL_filtered.vcf.gz


# selecting pass variants
zcat HCT116_DMSO_48h_SNP_filtered.vcf.gz | grep "#" > HCT116_DMSO_48h_SNP_filtered_PASS.vcf
zcat HCT116_DMSO_48h_INDEL_filtered.vcf.gz | grep "#" > HCT116_DMSO_48h_INDEL_filtered_PASS.vcf

zcat HCT116_DMSO_48h_SNP_filtered.vcf.gz | grep -v "#" | grep "PASS" >> HCT116_DMSO_48h_SNP_filtered_PASS.vcf
zcat HCT116_DMSO_48h_INDEL_filtered.vcf.gz | grep -v "#" | grep "PASS" >> HCT116_DMSO_48h_INDEL_filtered_PASS.vcf

# Merge
java -Xmx70g -jar ~/myPrograms/GenomeAnalysisTK.jar \
-T CombineVariants \
-R ~/resources/hg38/star/genome.fa \
--variant HCT116_DMSO_48h_SNP_filtered_PASS.vcf \
--variant HCT116_DMSO_48h_INDEL_filtered_PASS.vcf \
--genotypemergeoption UNIQUIFY \
-o HCT116_DMSO_48h_filtered_PASS.vcf

# Correcting fasta
java -Xmx70g -jar ~/myPrograms/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker \
    -R ~/resources/hg38/star/genome.fa  \
    -o /home/rtm/CY/RawData/P007_48_DMSO/HCT116.fasta \
--variant HCT116_DMSO_48h_filtered_PASS.vcf
