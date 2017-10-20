#!/bin/bash

#Bam files will be created using bwa mem
#Sample command: 
bwa mem -R \"@RG\tID:${ERSNO}\tSM:${SNAME}\" <GRCm38_68.fa|Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa>  $line.fq1.gz $line.fq2.gz | samtools view -bS - | samtools sort -o $line.bam  -O bam -T $line.temp

#Required if working with mouse exome sequencing data: bed files of all bait target regions separated by chromosome

#Sample progress of an analysis
filelist=''
while IFS='' read -r line || [[ -n "$line" ]]; do
  	n=$(echo $line | sed 's/.bam//g');
  	
  	#Run the samtools mpileup variant caller
  	samtools mpileup -f <GRCm38_68.fa|Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa> -g -t DP,DV -C50 -pm3 -F0.2 -d10000 $line | bcftools call -vm -f GQ > $n.vcf
  	
  	#If mouse run the scalpel variant caller chromosome by chromosome (this step can be paralleled with a cluster)
  	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr10_$n --bed ReferenceFiles/Bait_regions/chr10_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr11_$n --bed ReferenceFiles/Bait_regions/chr11_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr12_$n --bed ReferenceFiles/Bait_regions/chr12_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr13_$n --bed ReferenceFiles/Bait_regions/chr13_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr14_$n --bed ReferenceFiles/Bait_regions/chr14_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr15_$n --bed ReferenceFiles/Bait_regions/chr15_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr16_$n --bed ReferenceFiles/Bait_regions/chr16_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr17_$n --bed ReferenceFiles/Bait_regions/chr17_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr18_$n --bed ReferenceFiles/Bait_regions/chr18_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr19_$n --bed ReferenceFiles/Bait_regions/chr19_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr1_$n --bed ReferenceFiles/Bait_regions/chr1_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr2_$n --bed ReferenceFiles/Bait_regions/chr2_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr3_$n --bed ReferenceFiles/Bait_regions/chr3_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr4_$n --bed ReferenceFiles/Bait_regions/chr4_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr5_$n --bed ReferenceFiles/Bait_regions/chr5_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr6_$n --bed ReferenceFiles/Bait_regions/chr6_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr7_$n --bed ReferenceFiles/Bait_regions/chr7_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr8_$n --bed ReferenceFiles/Bait_regions/chr8_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir Chr9_$n --bed ReferenceFiles/Bait_regions/chr9_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	scalpel --somatic --normal Control.bam --tumor $line  --dir ChrX_$n --bed ReferenceFiles/Bait_regions/chrX_bait_region_mouse.bed --ref ReferenceGenomes/Mouse/GRCm38_68.fa
	#Concatenate the bam files
	#CONCATENATING FILES
	#Mutation files
	vcf-concat Chr1_$n/main/somatic.5x.indel.vcf Chr2_$n/main/somatic.5x.indel.vcf Chr3_$n/main/somatic.5x.indel.vcf Chr4_$n/main/somatic.5x.indel.vcf Chr5_$n/main/somatic.5x.indel.vcf Chr6_$n/main/somatic.5x.indel.vcf Chr7_$n/main/somatic.5x.indel.vcf Chr8_$n/main/somatic.5x.indel.vcf Chr9_$n/main/somatic.5x.indel.vcf Chr10_$n/main/somatic.5x.indel.vcf Chr11_$n/main/somatic.5x.indel.vcf Chr12_$n/main/somatic.5x.indel.vcf Chr13_$n/main/somatic.5x.indel.vcf Chr14_$n/main/somatic.5x.indel.vcf Chr15_$n/main/somatic.5x.indel.vcf Chr16_$n/main/somatic.5x.indel.vcf Chr17_$n/main/somatic.5x.indel.vcf Chr18_$n/main/somatic.5x.indel.vcf Chr19_$n/main/somatic.5x.indel.vcf ChrX_$n/main/somatic.5x.indel.vcf > $n.somatic.vcf
	#Common files
	vcf-concat Chr1_$n/main/common.5x.indel.vcf Chr2_$n/main/common.5x.indel.vcf Chr3_$n/main/common.5x.indel.vcf Chr4_$n/main/common.5x.indel.vcf Chr5_$n/main/common.5x.indel.vcf Chr6_$n/main/common.5x.indel.vcf Chr7_$n/main/common.5x.indel.vcf Chr8_$n/main/common.5x.indel.vcf Chr9_$n/main/common.5x.indel.vcf Chr10_$n/main/common.5x.indel.vcf Chr11_$n/main/common.5x.indel.vcf Chr12_$n/main/common.5x.indel.vcf Chr13_$n/main/common.5x.indel.vcf Chr14_$n/main/common.5x.indel.vcf Chr15_$n/main/common.5x.indel.vcf Chr16_$n/main/common.5x.indel.vcf Chr17_$n/main/common.5x.indel.vcf Chr18_$n/main/common.5x.indel.vcf Chr19_$n/main/common.5x.indel.vcf ChrX_$n/main/common.5x.indel.vcf > $n.common.vcf

	#Run the VEP on the samtools mpileup vcf files
	cat $n.vcf | variant_effect_predictor.pl --plugin Blosum62 --symbol --format vcf --cache --dir ref/ensembl/vep_cache --no_progress --quiet --offline --force_overwrite --species mus_musculus|saccharomyces_cerevisiae -o $n.vep.out
	vcf2consequences_vep -v $x -i $n.vep.out | bgzip -c > $n.cons.vcf.gz

	#Filter 
	#If mouse exon sequencing sample, filter non-exonic regions using the Bait_regions bed files
	#Scalpel somatic files
	cat $n.somatic.vcf | vcf-annotate -f gt-filter-lax.pl > $n.filt.vcf
	#Samtools mpileup result files
	cat $n.ex.vcf | vcf-annotate -H -f +/q=25/SnpGap=7/d=5 > $n.annotate.vcf
	cat $n.annotate.vcf |  vcf-annotate -f gt-filter-lax.pl > $n.gq.vcf
	
	#Normalising Indels
	bcftools norm -f <GRCm38_68.fa|Saccharomyces_cerevisiae.EF4.69.dna_sm.toplevel.fa>  $n.gq.vcf > $n.norm.vcf
	
	#Intersections of sample mpileup vcfs with Control mpileup vcfs
	cat $n.norm.vcf | vcf-sort > $n.sort.vcf
	bgzip -f $n.sort.vcf
	tabix -f -p vcf $n.sort.vcf.gz
	vcf-isec -f -a -c $x Control1.sort.vcf.gz Control2.sort.vcf.gz .. Controln.sort.vcf.gz > $n.isec.vcf

	#If mouse exons intersect samtools and scalpel results to retain only INDELs identified by both
	bedtools intersect -header -a $n.isec.vcf -b $n*somatic.filt.vcf > $n.merged.vcf
    cat $n.isec.vcf | grep "#" -v | grep "INDEL" -v >> $n.merged.vcf
	
	#Merge all files
	cat <$n.merged.vcf|$n.isec.vcf> | vcf-sort > $n.sort2.vcf 
	bgzip -f $n.sort2.vcf 
	bgzip -f $n.sort2.vcf 
	tabix -p vcf $n.sort2.vcf.gz	
	
	#create a list of all final vcf files created
	b='.sort2.vcf.gz '
	c=$n$b
	filelist+=$c

done < "$1" #provide a list of bam files from the command line (with one bam file per line)

#Continue merge
vcf-merge $filelist  > Experiment_merge.vcf

#Create the list of potential suppressor mutations
perl vcf_to_gene_list.pl -i Experiment_merge.vcf > results.tsv

