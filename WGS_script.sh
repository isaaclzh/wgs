# Germline Short Variant Calling

#!/usr/bin/env sh

# Working folder
FOLDER_NAME="ashkenazim_giab_results" #This will be your working directory which contains all the output files
FOLDER="ashkenazim_samples" #This is the folder which contains your fastq files (Home Directory)
gatk=/home/zhenhanisaac.lin/gatk-4.3.0.0/gatk #The is the location which contains your gatk executable
i=0 # Looping variable

# Load modules
module load bwa
module load samtools
module load bcftools
module load bamtools
module load python
module load trimgalore
module load r
module load vep

# Error catch
set -euo pipefail

# Create directories
mkdir -p ~/$FOLDER_NAME && cd $_
mkdir -p bam
mkdir -p vcf/gvcf
mkdir -p reports/plot
mkdir -p reports/duplicate_metrics

cd ~/$FOLDER

# BWA-MEM alignment
for f in `ls -1 *_1.fq.gz | sed 's/_1.fq.gz//'`; do
	echo "Performing BWA MEM on ${f}_1.fq.gz ${f}_2.fq.gz"
	bwa mem -M -t 20 \
		-R "@RG\tID:${f}\tSM:${f}\tPL:Illumina\tLB:${f}" \
		~/reference_genome/GRCh38_Verily_v1.genome.fa ${f}_1.fq.gz ${f}_2.fq.gz \
		| samtools view -Sb - \
		| samtools sort - -@ 20 -m 4G -o ${f}.bwa.bam
done

# Moving files
mv `ls *.bam` ~/$FOLDER_NAME/bam

cd ~/$FOLDER_NAME/bam

# Create a cohort file for haplotypecaller
:> cohort.txt

# Remove duplicates
for f in `ls -1 *.bam | sed 's/.bam'//`; do
	echo "Marking duplicates on ${f}.bam"	
	$gatk --java-options "-Xmx20G -XX:ParallelGCThreads=8" MarkDuplicates \
		-I ${f}.bam \
        	-O ${f}.sortdup.bam \
        	-M ${f}.sortdup.bam.duplicate_metrics \
		--VALIDATION_STRINGENCY SILENT \
        	--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        	--ASSUME_SORT_ORDER coordinate \
        	--CLEAR_DT false \
        	--REMOVE_DUPLICATES true \
        	--CREATE_INDEX false \
		--READ_NAME_REGEX null

        mv ${f}.sortdup.bam.duplicate_metrics ../reports/duplicate_metrics
done

# Base recalibrator
for f in `ls -1 *.bam | sed 's/.bam'//`; do
	echo "Running base recalibrator on ${f}.bam"
	$gatk --java-options "-Xmx4G" BaseRecalibrator \
		-I ${f}.bam \
		-O ${f}.recal_data.table \
		-R ~/reference_genome/GRCh38_Verily_v1.genome.fa \
		-OQ true \
		--known-sites ~/reference_genome/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
		--known-sites ~/reference_genome/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
		--known-sites ~/reference_genome/Homo_sapiens_assembly38.known_indels.vcf.gz

	$gatk --java-options "-Xmx4G" ApplyBQSR \
		-I ${f}.bam \
		-O ${f}.bqsr.bam \
		-R ~/reference_genome/GRCh38_Verily_v1.genome.fa \
		-bqsr ${f}.recal_data.table \
		-OQ true \
		-OBI false \
		--add-output-sam-program-record true

	samtools index ${f}.bqsr.bam
done

# HaplotypeCaller
for f in `ls -1 *.bwa.sortdup.bqsr.bam | sed 's/.bwa.sortdup.bqsr.bam'//`; do
	echo "Running GATK HaplotypeCaller on ${f}.bam"
	$gatk --java-options "-Xmx16G" HaplotypeCaller \
		--native-pair-hmm-threads 4 \
		-I ${f}.bwa.sortdup.bqsr.bam \
		-R ~/reference_genome/GRCh38_Verily_v1.genome.fa \
		-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 \
                -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY \
		-O ${f}.g.vcf.gz \
		-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
		--pcr-indel-model NONE \
		-ERC GVCF
	
	# Writing samples
	i=`expr $i + 1`
	echo -e "sample $i\t${f}.g.vcf.gz" >> cohort.txt
done

# Moving files
mv cohort.txt ../vcf
mv `ls *.vcf.*` ../vcf/gvcf

# Delete all intermediate bam files
find . -type f ! -name '.bam.' -delete

cd ../vcf/gvcf

# GenomicsDBImport
$gatk --java-options "-Xmx4G" GenomicsDBImport \
	--genomicsdb-workspace-path ../database \
	--batch-size 0 \
	-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 \
        -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY \
	--sample-name-map ../cohort.txt \
	--reader-threads 5

# GenotypeGVCFs
$gatk --java-options "-Xmx5G" GenotypeGVCFs \
	-R ~/reference_genome/GRCh38_Verily_v1.genome.fa \
	-O ../raw_variants.vcf.gz \
	-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
	-V gendb://../database \
	-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 \
        -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY \
	--only-output-calls-starting-in-intervals true \
	--use-new-qual-calculator true

cd ..

# VariantRecalibrator (SNP)
$gatk --java-options "-Xmx16G" VariantRecalibrator \
	-R ~/reference_genome/GRCh38_Verily_v1.genome.fa \
	-V raw_variants.vcf.gz \
	-O recalibrate_snp.recal \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	--resource:hapmap,known=false,training=true,truth=true,prior=15 ~/reference_genome/hapmap_3.3.hg38.vcf.gz \
	--resource:omni,known=false,training=true,truth=false,prior=12 ~/reference_genome/1000G_omni2.5.hg38.vcf.gz \
	--resource:1000G,known=false,training=true,truth=false,prior=10 ~/reference_genome/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	--resource:dbsnp,known=true,training=false,truth=false,prior=7 ~/reference_genome/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
	-mode SNP \
	--max-gaussians 6 \
	--tranches-file recalibrate_snp.tranches \
	--rscript-file plots_snp.R

# ApplyVQSR (SNP)
$gatk --java-options "-Xmx16G" ApplyVQSR \
	-R ~/reference_genome/GRCh38_Verily_v1.genome.fa \
	-V raw_variants.vcf.gz \
	-O recalibrated_snps_raw_indels.vcf.gz \
	-mode SNP \
	-ts-filter-level 99.7 \
	--tranches-file recalibrate_snp.tranches \
	--recal-file recalibrate_snp.recal \
	--create-output-variant-index true

# VariantRecalibrator (INDEL)
$gatk --java-options "-Xmx16G" VariantRecalibrator \
	-R ~/reference_genome/GRCh38_Verily_v1.genome.fa \
	-V recalibrated_snps_raw_indels.vcf.gz \
	-O recalibrate_indel.recal \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	--resource:mills,known=false,training=true,truth=true,prior=12 ~/reference_genome/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	--resource:axiomPoly,known=false,training=true,truth=false,prior=10 ~/reference_genome/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2 ~/reference_genome/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
	-mode INDEL \
	--max-gaussians 4 \
	--tranches-file recalibrate_indel.tranches \
	--rscript-file plots_indel.R

# ApplyVQSR (INDEL)
$gatk --java-options "-Xmx5G" ApplyVQSR \
	-R ~/reference_genome/GRCh38_Verily_v1.genome.fa \
	-V recalibrated_snps_raw_indels.vcf.gz \
	-O recalibrated_variants.vcf.gz \
	--mode INDEL \
	-ts-filter-level 99.7 \
	--tranches-file recalibrate_indel.tranches \
	--recal-file recalibrate_indel.recal \
	--create-output-variant-index true

chmod +x *.vcf.gz

# Moving Files
mv `ls plots*` ../reports/plot

<<C
# Step 5: Variant Annotation
echo "Running Variant Effect Predictor"
vep -i recalibrated_variants.vcf.gz \
	-o ${FOLDER_NAME}_annotation.txt \ 
	--refseq \
	--use_given_ref \
	--hgvs \
	--fasta /home/cbi/biodata/gencode/GRCh38/GRCh38.p13.genome.fa
C
