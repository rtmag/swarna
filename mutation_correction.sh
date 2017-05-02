
mkdir mutation_correction
cd mutation_correction

# 2.1.1) Variant Calling

# BWA installation
git clone https://github.com/lh3/bwa.git
cd bwa
ln -s /home/rtm/myPrograms/bwa-0.7.15/bwa /home/rtm/bin/bwa

# BWA generating genome index
bwa index -a bwtsw fly.fa

# BWA mapping
bwa mem -M -t 16 fly.fa wgs_1.fq.gz wgs_2.fq.gz > wgs.sam
# -M ensures compatibility with picard
# -t number of threads

# SAM TO BAM + SORT
samtools view -Sb wgs.sam |samtools sort - wgs

# Fasta index
samtools faidx genome.fa
# Creating genome dir
java -jar picard.jar CreateSequenceDictionary R= fly.fa O= genome.dict
# Add RG; CHECK FOR RG IN BAM HEADER IF IS THERE NO NEED FOR THIS STEP
java -jar picard.jar AddOrReplaceReadGroups INPUT=wgs.bam OUTPUT=wgs_addRG.bam RGID=HNHCCCCXX RGLB= Merged RGPL=illumina RGPU=HNHCCCCXX RGSM=sample1 
# Index
samtools index wgs_addRG.bam

# ReMapping
java -Xmx50g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R fly.fa \
-I wgs_addRG.bam -o lane.intervals 
--known fly_indels.vcf

# INDEL ReMapping
java -Xmx60g -jar GenomeAnalysisTK.jar -T IndelRealigner -nct 14 -R fly.fa \
-I wgs_addRG.bam -targetIntervals lane.intervals \
-known fly_indels.vcf -o wgs_addRG_realigned.bam

#Base recallibration 
java -Xmx60g -jar GenomeAnalysisTK.jar -T BaseRecalibrator -nct 14 -R fly.fa -knownSites fly_indels.vcf -I wgs_addRG_realigned.bam -o wgs_recal.table
java -Xmx60g -jar GenomeAnalysisTK.jar -T PrintReads -nct 14 -R fly.fa -I wgs_addRG_realigned.bam --BQSR wgs_recal.table -o  wgs_addRG_realigned_recalibrated.bam

# Remove duplicates
java -Xmx60g -jar picard.jar MarkDuplicates  VALIDATION_STRINGENCY=LENIENT M=MarkDup_metrics.txt INPUT=wgs_addRG_realigned_recalibrated.bam OUTPUT=wgs_addRG_realigned_recalibrated_rmdup.bam

#INDEX
samtools index wgs_addRG_realigned_recalibrated_rmdup.bam

# Variant calling pileup IM USING samtools 1.3.1
samtools mpileup -ugf fly.fa wgs_addRG_realigned_recalibrated_rmdup.bam | bcftools call -vmO z -o wgs.vcf.gz

# VCF index
tabix -p vcf wgs.vcf.gz

# Extract SNP and INDELS
java -jar GenomeAnalysisTK.jar -T SelectVariants \
-R fly.fa \
-V wgs.vcf.gz \
-selectType SNP \
-o wgs_SNP.vcf.gz 

# Filtering

java -jar GenomeAnalysisTK.jar \
-T VariantFiltration \
-R fly.fa \
-V wgs_SNP.vcf.gz \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filterName "snp_filter" \
-o wgs_SNP_filtered.vcf.gz 


# selecting pass variants
zcat wgs_SNP_filtered.vcf.gz | grep -v "#" | grep "PASS" >> wgs_SNP_filtered_PASS.vcf

# Correcting fasta
java -Xmx70g -jar GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker \
    -R fly.fa  \
    -o fly_corrected.fasta \
--variant wgs_filtered_PASS.vcf
