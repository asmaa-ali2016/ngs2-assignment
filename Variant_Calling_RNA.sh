###########################################
#  - Variants Calling in RNAseq           #
#  - NGS2 Course - Assignment             #
#  - Bash Script                          #
#  - April 27,2019                        #
#  - Copyright: Asmaa Ali                 #
#  - Nile University                      #
###########################################

#!/bin/bash


# Install Packages
conda install -c bioconda star 
conda install -c bioconda samtools 
conda install -c bioconda picard 
conda install -c bioconda gatk4 
conda install -c bioconda tabix
conda install -c bioconda rtg-tools
sudo apt install picard-tools

cd ~/NU_Bioinformatics_Diploma/NGS2/Assignment

# prepare paths
STAR=/home/dna/STAR/bin/Linux_x86_64/STAR 
GATK=/home/dna/anaconda3/share/gatk4-4.1.1.0-0/gatk-package-4.1.1.0-local.jar
GenomeDir=/home/dna/NU_Bioinformatics_Diploma/NGS2/Assignment/index
GenomeFasta=/home/dna/NU_Bioinformatics_Diploma/NGS2/Assignment/index/Homo_sapiens.GRCh38.dna.chromosome.15.fa
Reads="/home/dna/NU_Bioinformatics_Diploma/NGS2/Assignment/fqData/*.fastq"
CommonPars="--runThreadN 8 --outSAMattributes All"

# generate indexed genome for 1st pass alignment
mkdir STAR_index
cd STAR_index
$STAR --genomeDir $GenomeDir --runMode genomeGenerate --genomeFastaFiles $GenomeFasta \
--sjdbGTFfile /home/dna/NU_Bioinformatics_Diploma/NGS2/Assignment/index/Homo_sapiens.GRCh38.dna.chromosome.15.gtf --sjdbOverhang 100 --runThreadN 8

# run 1st pass
mkdir Pass1
cd Pass1
for f in $Reads; do \
$STAR $CommonPars --genomeDir $GenomeDir --readFilesIn $f --outFileNamePrefix ${f}_pass1 ; done
cd .. && mv fqData/*pass1* Pass1/

# make splice junctions database file out of SJ.out.tab, 
# filter out non-canonical junctions 
mkdir GenomeForPass2
cd GenomeForPass2
awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} \
{if($5 > 0){print $1,$2,$3,strChar[$4]}}' \
../Pass1/*SJ.out.tab > SJ.out.tab.Pass1.sjdb

# generate genome with junctions from the 1st pass
$STAR --genomeDir ./ --runMode genomeGenerate --genomeFastaFiles $GenomeFasta \
--sjdbFileChrStartEnd SJ.out.tab.Pass1.sjdb --sjdbOverhang 100 --runThreadN 8
cd ..

# run 2nd pass with the new genome (for novel junction discovery)
mkdir Pass2
cd Pass2
for f in $Reads; do \ 
$STAR $CommonPars --genomeDir ../GenomeForPass2 --readFilesIn $f \
--outFileNamePrefix ${f}_pass2 ; done
cd .. && mv fqData/*pass2* Pass2/


# generate & sort BAM file
for samfile in Pass2/*.sam;do
  sample=${samfile%.sam}
  samtools view -hbo $sample.bam $samfile
  samtools sort $sample.bam -o $sample.sorted.bam
done

# mapping QC
for bamFile in Pass2/*.sorted.bam;do
  output=${bamFile%.sorted.bam}
  samtools depth $bamFile | awk '{{sum+=$3}} END {{print "Average = ",sum/NR}}' > $output.cov
  samtools flagstat $bamFile > $output.stat
done


# Picard Markduplicates
picard_path=/home/dna/anaconda3/share/picard-2.19.0-0
for sample in Pass2/*.sorted.bam;do
  name=${sample%.sorted.bam}
  java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates INPUT=$sample OUTPUT=$name.dedup.bam METRICS_FILE=$name.metrics.txt;
done


# index bam files   
for b in Pass2/*.dedup.bam ; do
samtools index $b ; done

# Reference Sequence Dictionary file (AAA.dict)
picard-tools CreateSequenceDictionary R=$GenomeFasta O=$GenomeFasta.dict
samtools faidx $GenomeFasta
#mv $GenomeFasta.fai index/Homo_sapiens.GRCh38.dna.chromosome.15.fai
#mv $GenomeFasta.dict index/Homo_sapiens.GRCh38.dna.chromosome.15.dict


# Split'N'Trim and reassign mapping qualities
for d in Pass2/*.dedup.bam ; do \
java -jar $GATK SplitNCigarReads -R $GenomeFasta -I $d -O $d.split.bam ; done


# Recalibrate Bases BQSR
for sample in Pass2/*.split.bam;do
  gatk --java-options "-Xmx2G" BaseRecalibrator \
-R $GenomeFasta -I $sample --known-sites Homo_sapiens.vcf \
-O $sample.report

  gatk --java-options "-Xmx2G" ApplyBQSR \
-R $GenomeFasta -I $sample -bqsr $sample.report \
-O $sample.bqsr.bam --add-output-sam-program-record --emit-original-quals
done


for sample in Pass2/*.bqsr.bam;do
samtools view $sample | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq > ${sample}_samplename.txt
done


## assess genotype likelihood per-sample
for sample in Pass2/*.bqsr.bam;do
  name=${sample%.bqsr.bam}

  gatk --java-options "-Xmx2G" HaplotypeCaller \
  -I $sample -O $name.gvcf\
  --emit-ref-confidence GVCF \
  --pcr-indel-model NONE \
  -R $GenomeFasta
done 


## combine samples
gatk --java-options "-Xmx2G" CombineGVCFs \
-R $GenomeFasta \
-V SRR8797509_1.part_001.part_001_pass2Aligned.out.dedup.split.bqsr.gvcf \
-V SRR8797509_2.part_001.part_001_pass2Aligned.out.dedup.split.bqsr.gvcf \
-O raw_variants.gvcf


## Joint Genotyping
gatk --java-options "-Xmx60G" GenotypeGVCFs \
-R $GenomeFasta \
-V raw_variants.gvcf \
--max-alternate-alleles 6 \
-O raw_variants.vcf


## annotated output
gatk --java-options "-Xmx60G" GenotypeGVCFs \
-R $GenomeFasta \
-V raw_variants.gvcf \
--max-alternate-alleles 6 \
--dbsnp canis_fam_chr5.vcf \
-O raw_variants_ann.vcf


## check how many variant got annotated
grep -v "^#" raw_variants_ann.vcf | awk '{print $3}' | grep "^rs" | wc -l

# VCF statitics
bgzip -c raw_variants_ann.vcf > raw_variants_ann.vcf.gz
tabix -p vcf raw_variants_ann.vcf.gz
rtg vcfstats raw_variants_ann.vcf.gz > stats.txt


# Split SNPs and indels
gatk --java-options "-Xmx2G" SelectVariants \
-R dog_chr5.fa \
-V raw_variants_ann.vcf \
--select-type-to-include SNP \
-O raw_variants_ann_SNP.vcf

gatk --java-options "-Xmx2G" SelectVariants \
-R dog_chr5.fa \
-V raw_variants_ann.vcf \
--select-type-to-include INDEL \
-O raw_variants_ann_INDEL.vcf


# Assess the different filters in both known and novel
for var in "SNP" "INDEL";do
 input="raw_variants_ann_"$var".vcf"
 for filter in "QD" "MQ" "MQRankSum" "FS" "SOR" "ReadPosRankSum" "AN" "DP" "InbreedingCoeff";do
  filterValues=$var.$filter
  awk -v k="$filter=" '!/#/{n=split($8,a,";"); for(i=1;i<=n;i++) if(a[i]~"^"k) {sub(k,$3" ",a[i]); print a[i]}}' $input > $filterValues
  grep -v "^\." $filterValues > known.$var.$filter
  grep "^\." $filterValues > novel.$var.$filter
done; done
mkdir filters && cd filters
mv ../{*.SNP.*,SNP.*,*.INDEL.*,INDEL.*} .


# Figures
wget https://raw.githubusercontent.com/dib-lab/dogSeq/master/scripts/densityCurves.R
sudo Rscript -e "install.packages('ggplot2', contriburl=contrib.url('http://cran.r-project.org/'))"
for f in SNP.* INDEL.*;do
  Rscript densityCurves.R "$f"
done


# Calc the DP threathols
cat SNP.DP INDEL.DP | awk '{sum+= $2; sumsq+= ($2)^2} END { print sum/NR, sqrt((sumsq-sum^2/NR)/NR), sum/NR + 5*sqrt((sumsq-sum^2/NR)/NR) }' 


# SNP Variant filteration
for v in Pass2/*.bam.vcf; do 
java -jar $GATK \
VariantFiltration \
-R $GenomeFasta \
-V $v \
--filter-name "snpQD" \
--filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
--filter-name "snpMQ" \
--filter-expression "vc.hasAttribute('MQ') && MQ < 40.0" \
--filter-name "snpMQRankSum" \
--filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
--filter-name "snpFS" \
--filter-expression "vc.hasAttribute('FS') && FS > 60.0" \
--filter-name "snpSOR" \
--filter-expression "vc.hasAttribute('SOR') && SOR > 4.0" \
--filter-name "snpReadPosRankSum" \
--filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
--filter-name "snpDP" \
--filter-expression "vc.hasAttribute('DP') && DP > 3105" \
-O $v.snp.filtered.vcf; done


# INDEL Variant filteration
for v in Pass2/*.bam.vcf; do 
java -jar $GATK \
VariantFiltration \
-R $GenomeFasta \
-V $v \
--filter-name "indelQD" \
--filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
--filter-name "indelFS" \
--filter-expression "vc.hasAttribute('FS') && FS > 200.0" \
--filter-name "indelSOR" \
--filter-expression "vc.hasAttribute('SOR') && SOR > 10.0" \
--filter-name "indelReadPosRankSum" \
--filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \
--filter-name "indelInbreedingCoeff" \
--filter-expression "vc.hasAttribute('InbreedingCoeff') && InbreedingCoeff < -0.8" \
--filter-name "indelDP" \
--filter-expression "vc.hasAttribute('DP') && DP > 3105" \
-O $v.indel.filtered.vcf; done


# vcf files can be uploaded on IGV
