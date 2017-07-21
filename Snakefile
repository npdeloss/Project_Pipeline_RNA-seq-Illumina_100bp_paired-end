import os
import glob

bin_home_dir = "/home/npdeloss/bin"

samtools = bin_home_dir + "/samtools"
STAR = bin_home_dir + "/STAR"
fastqc = bin_home_dir + "/fastqc"
rsem_prepare_reference = bin_home_dir + "/rsem-prepare-reference"
rsem_calculate_expression = bin_home_dir + "/rsem-calculate-expression"
gtf2bed = bin_home_dir + "/gtf2bed"
rseqc_bin = bin_home_dir

processors_per_node = 16
email_address = "npdeloss@ucsd.edu"
cluster_computing_group = "jogleeson-group"

human_sample_names = list(set([path.split('/')[-1].split('-sample_concat-')[0] for path in glob.glob('Raw_fastq_gz/H_sapiens/*-sample_concat-*_R1_*.fastq.gz')]))

rule all:
	input:
		expand("QC_fastqc/H_sapiens/{sample}_R1", sample = human_sample_names),
		expand("Alignment_STAR/H_sapiens/GRCh38/{sample}/Aligned.out.bam", sample = human_sample_names),
		expand("QC_RSeQC/H_sapiens/GRCh38/{sample}", sample = human_sample_names),
		expand("Quantification_RSEM/H_sapiens/GRCh38/{sample}/Quant.genes.results", sample = human_sample_names)
		
rule download_reference_GRCh38:
	output:
		reference_dir = "Reference/H_sapiens/GRCh38/",
		reference_fasta = "Reference/H_sapiens/GRCh38/GRCh38.primary_assembly.genome.fa",
		reference_gtf = "Reference/H_sapiens/GRCh38/gencode.v26.annotation.gtf"
	params:
		cluster = " -N download_reference_GRCh38 "
		"-o  Reference/H_sapiens/GRCh38_job.log "
		"-e  Reference/H_sapiens/GRCh38_job.err "
		"-q hotel "
		"-l nodes=1:ppn=1 "
		"-l walltime=1:00:00 "
		"-V " +
		"-M " + email_address + " "
		"-m abe "
		"-A " + cluster_computing_group
	shell:
		"mkdir -p Reference/H_sapiens/GRCh38/ ;"
		"cd Reference/H_sapiens/GRCh38/ ;"
		"wget -nc ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/GRCh38.primary_assembly.genome.fa.gz ;"
		"zcat GRCh38.primary_assembly.genome.fa.gz > GRCh38.primary_assembly.genome.fa ;"
		"wget -nc ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz ;"
		"zcat gencode.v26.annotation.gtf.gz > gencode.v26.annotation.gtf ;"

rule STAR_reference_GRCh38:
	input:
		reference_fasta = "Reference/H_sapiens/GRCh38/GRCh38.primary_assembly.genome.fa",
		reference_gtf = "Reference/H_sapiens/GRCh38/gencode.v26.annotation.gtf"
	output:
		STAR_reference_dir = "Reference/H_sapiens/GRCh38/STAR"
	params:
		cluster = " -N STAR_reference_GRCh38 "
		"-o Reference/H_sapiens/GRCh38/STAR_job.log "
		"-e Reference/H_sapiens/GRCh38/STAR_job.err "
		"-q hotel "
		"-l nodes=1:ppn=" + str(processors_per_node) + " "
		"-l walltime=24:00:00 "
		"-V "
		"-M " + email_address + " "
		"-m abe "
		"-A " + cluster_computing_group
	shell:
		"mkdir -p {output.STAR_reference_dir} ;"
		"cd {output.STAR_reference_dir} ;"
		"{STAR} "
		"--runThreadN {processors_per_node} "
		"--runMode genomeGenerate "
		"--genomeDir ./ "
		"--genomeFastaFiles ../GRCh38.primary_assembly.genome.fa "
		"--sjdbGTFfile ../gencode.v26.annotation.gtf "
		"--sjdbOverhang 99"

rule RSEM_reference_GRCh38:
	input:
		reference_fasta = "Reference/H_sapiens/GRCh38/GRCh38.primary_assembly.genome.fa",
		reference_gtf = "Reference/H_sapiens/GRCh38/gencode.v26.annotation.gtf"
	output:
		RSEM_reference_dir = "Reference/H_sapiens/GRCh38/RSEM"
	params:
		cluster = " -N RSEM_reference_GRCh38 "
		"-o Reference/H_sapiens/GRCh38/RSEM_job.log "
		"-e Reference/H_sapiens/GRCh38/RSEM_job.err "
		"-q hotel "
		"-l nodes=1:ppn=1" + " "
		"-l walltime=24:00:00 "
		"-V "
		"-M " + email_address + " "
		"-m abe "
		"-A " + cluster_computing_group
	shell:
		"mkdir -p {output.RSEM_reference_dir} ;"
		"cd {output.RSEM_reference_dir} ;"
		"{rsem_prepare_reference}  --gtf $(ls ../*.gtf) $(ls ../*.fa) reference"

rule RSeQC_reference_GRCh38:
	input:
		reference_gtf = "Reference/H_sapiens/GRCh38/gencode.v26.annotation.gtf"
	output:
		RSeQC_reference_dir = "Reference/H_sapiens/GRCh38/RSeQC",
		RSeQC_reference_bed = "Reference/H_sapiens/GRCh38/RSeQC/gene_model.bed"
	params:
		cluster = " -N RSeQC_reference_GRCh38 "
		"-o Reference/H_sapiens/GRCh38/RSeQC_job.log "
		"-e Reference/H_sapiens/GRCh38/RSeQC_job.err "
		"-q hotel "
		"-l nodes=1:ppn=1" + " "
		"-l walltime=24:00:00 "
		"-V "
		"-M " + email_address + " "
		"-m abe "
		"-A " + cluster_computing_group
	shell:
		"mkdir -p {output.RSeQC_reference_dir} ;"
		"{gtf2bed} {input.reference_gtf} > {output.RSeQC_reference_bed} ;"

rule Concatenate_fastq_gz:
	input:
		Raw_fastq_gz_R1 = lambda wildcards: glob.glob(expand("Raw_fastq_gz/{organism}/{sample}-sample_concat-*_R1_*.fastq.gz",organism=wildcards.organism,sample=wildcards.sample)[0]),
		Raw_fastq_gz_R2 = lambda wildcards: glob.glob(expand("Raw_fastq_gz/{organism}/{sample}-sample_concat-*_R2_*.fastq.gz",organism=wildcards.organism,sample=wildcards.sample)[0])
	output:
		Concatenate_fastq_gz_R1 = "Concatenate_fastq_gz/{organism}/{sample}_R1.fastq.gz",
		Concatenate_fastq_gz_R2 = "Concatenate_fastq_gz/{organism}/{sample}_R2.fastq.gz"
	params:
		cluster = lambda wildcards: expand(" -N Concatenate_fastq_gz_{organism}_{sample} " +
		"-o Concatenate_fastq_gz/{organism}/{sample}_job.log " +
		"-e Concatenate_fastq_gz/{organism}/{sample}_job.err " +
		"-q hotel " +
		"-l nodes=1:ppn=2" + " " +
		"-l walltime=24:00:00 " +
		"-V " +
		"-M " + email_address + " " +
		"-m abe " +
		"-A " + cluster_computing_group + "",
		organism=wildcards.organism, 
		sample=wildcards.sample)[0]
	shell:
		"""mkdir -p Concatenate_fastq_gz/{wildcards.organism}
		zcat Raw_fastq_gz/{wildcards.organism}/{wildcards.sample}-sample_concat-*_R1_*.fastq.gz | gzip > Concatenate_fastq_gz/{wildcards.organism}/{wildcards.sample}_R1.fastq.gz &
		zcat Raw_fastq_gz/{wildcards.organism}/{wildcards.sample}-sample_concat-*_R2_*.fastq.gz | gzip > Concatenate_fastq_gz/{wildcards.organism}/{wildcards.sample}_R2.fastq.gz &
		wait"""

rule Alignment_STAR:
	input:
		STAR_reference_dir = "Reference/{organism}/{reference}/STAR",
		Concatenate_fastq_gz_R1 = "Concatenate_fastq_gz/{organism}/{sample}_R1.fastq.gz",
		Concatenate_fastq_gz_R2 = "Concatenate_fastq_gz/{organism}/{sample}_R2.fastq.gz"
	output:
		unsorted_BAM = "Alignment_STAR/{organism}/{reference}/{sample}/Aligned.out.bam",
		sorted_BAM = "Alignment_STAR/{organism}/{reference}/{sample}/Aligned.sortedByCoord.out.bam",
		transcriptome_BAM = "Alignment_STAR/{organism}/{reference}/{sample}/Aligned.toTranscriptome.out.bam"
	params:
		cluster = lambda wildcards: expand(" -N Alignment_STAR_{organism}_{reference}_{sample} " +
		"-o Alignment_STAR/{organism}/{reference}/{sample}_job.log " +
		"-e Alignment_STAR/{organism}/{reference}/{sample}_job.err " +
		"-q hotel " +
		"-l nodes=1:ppn=" + str(processors_per_node) + " " +
		"-l walltime=24:00:00 " +
		"-V " +
		"-M " + email_address + " " +
		"-m abe " +
		"-A " + cluster_computing_group + "",
		organism=wildcards.organism,
		reference=wildcards.reference,
		sample=wildcards.sample)[0]
	shell:
		"mkdir -p Alignment_STAR/{wildcards.organism}/{wildcards.reference}/{wildcards.sample} ;"
		"cd Alignment_STAR/{wildcards.organism}/{wildcards.reference}/{wildcards.sample} ;"
		"{STAR} "
		"--runThreadN {processors_per_node} "
		"--genomeDir ../../../../{input.STAR_reference_dir} "
		"--twopassMode Basic "
		"--sjdbGTFfile $(ls ../../../../{input.STAR_reference_dir}/../*.gtf) "
		"--readFilesIn ../../../../{input.Concatenate_fastq_gz_R1} ../../../../{input.Concatenate_fastq_gz_R2} "
		"--sjdbOverhang 99 "
		"--readFilesCommand zcat "
		"--outBAMcompression 9 "
		"--outSAMtype BAM Unsorted SortedByCoordinate "
		"--quantMode TranscriptomeSAM GeneCounts "
		"--outSAMattributes NH HI AS nM NM MD jM jI XS ch "
		"--outSAMunmapped Within "
		"--outReadsUnmapped Fastx "
		"--outFilterType BySJout "
		"--outFilterMultimapNmax 20 "
		"--outFilterMismatchNmax 999 "
		"--outFilterMismatchNoverLmax 0.04 "
		"--alignIntronMin 20 "
		"--alignIntronMax 1000000 "
		"--alignMatesGapMax 1000000 "
		"--alignSJoverhangMin 8 "
		"--alignSJDBoverhangMin 1 ;"
		"wait ;"
		"sleep 1 ;"
		"{samtools} index Aligned.sortedByCoord.out.bam ;"

rule Quantification_RSEM:
	input:
		RSEM_reference_dir = "Reference/{organism}/{reference}/RSEM",
		transcriptome_BAM = "Alignment_STAR/{organism}/{reference}/{sample}/Aligned.toTranscriptome.out.bam"
	output:
		isoform_quantifications = "Quantification_RSEM/{organism}/{reference}/{sample}/Quant.isoforms.results",
		gene_quantifications = "Quantification_RSEM/{organism}/{reference}/{sample}/Quant.genes.results"
	params:
		cluster = lambda wildcards: expand(" -N Quantification_RSEM_{organism}_{reference}_{sample} " +
		"-o Quantification_RSEM/{organism}/{reference}/{sample}_job.log " +
		"-e Quantification_RSEM/{organism}/{reference}/{sample}_job.err " +
		"-q hotel " +
		"-l nodes=1:ppn=" + str(processors_per_node) + " " +
		"-l walltime=24:00:00 " +
		"-V " +
		"-M " + email_address + " " +
		"-m abe " +
		"-A " + cluster_computing_group + "",
		organism=wildcards.organism,
		reference=wildcards.reference,
		sample=wildcards.sample)[0]
	shell:
		"mkdir -p Quantification_RSEM/{wildcards.organism}/{wildcards.reference}/{wildcards.sample} ;"
		"cd Quantification_RSEM/{wildcards.organism}/{wildcards.reference}/{wildcards.sample} ;"
		"{rsem_calculate_expression} "
		"--bam "
		"--estimate-rspd "
		"--no-bam-output "
		"--seed 12345 "
		"-p {processors_per_node} "
		"--paired-end --forward-prob 0.0 "
		"../../../../{input.transcriptome_BAM} "
		"../../../../{input.RSEM_reference_dir}/reference "
		"Quant "

rule QC_fastqc:
	input:
		Concatenate_fastq_gz_R1 = "Concatenate_fastq_gz/{organism}/{sample}_R1.fastq.gz",
		Concatenate_fastq_gz_R2 = "Concatenate_fastq_gz/{organism}/{sample}_R2.fastq.gz"
	output:
		QC_fastqc_R1 = "QC_fastqc/{organism}/{sample}_R1",
		QC_fastqc_R2 = "QC_fastqc/{organism}/{sample}_R2"
	params:
		cluster = lambda wildcards: expand(" -N QC_fastqc_{organism}_{sample} " +
		"-o QC_fastqc/{organism}/{sample}_job.log " +
		"-e QC_fastqc/{organism}/{sample}_job.err " +
		"-q hotel " +
		"-l nodes=1:ppn=1 " +
		"-l walltime=24:00:00 " +
		"-V " +
		"-M " + email_address + " " +
		"-m abe " +
		"-A " + cluster_computing_group + "",
		organism=wildcards.organism,
		sample=wildcards.sample)[0]
	shell:
		"mkdir -p QC_fastqc/{wildcards.organism} ;"
		"mkdir -p {output.QC_fastqc_R1}/tmp ;"
		"mkdir -p {output.QC_fastqc_R2}/tmp ;"
		"{fastqc} -o {output.QC_fastqc_R1} -d {output.QC_fastqc_R1}/tmp {input.Concatenate_fastq_gz_R1} ;"
		"{fastqc} -o {output.QC_fastqc_R2} -d {output.QC_fastqc_R2}/tmp {input.Concatenate_fastq_gz_R2} ;"
		"rm -rf {output.QC_fastqc_R1}/tmp ;"
		"rm -rf {output.QC_fastqc_R2}/tmp ;"

rule QC_RSeQC:
	input:
		unsorted_BAM = "Alignment_STAR/{organism}/{reference}/{sample}/Aligned.out.bam",
		sorted_BAM = "Alignment_STAR/{organism}/{reference}/{sample}/Aligned.sortedByCoord.out.bam",
		RSeQC_reference_bed = "Reference/{organism}/{reference}/RSeQC/gene_model.bed",
		reference_dir = "Reference/{organism}/{reference}/"
	output:
		"QC_RSeQC/{organism}/{reference}/{sample}"
	params:
		cluster = lambda wildcards: expand(" -N QC_RSeQC_{organism}_{sample} " +
		"-o QC_RSeQC/{organism}/{reference}/{sample}_job.log " +
		"-e QC_RSeQC/{organism}/{reference}/{sample}_job.err " +
		"-q hotel " +
		"-l nodes=1:ppn=16 " +
		"-l walltime=24:00:00 " +
		"-V " +
		"-M " + email_address + " " +
		"-m abe " +
		"-A " + cluster_computing_group + "",
		organism=wildcards.organism,
		reference=wildcards.reference,
		sample=wildcards.sample)[0]
	shell:
		"""mkdir -p QC_RSeQC/{wildcards.organism}/{wildcards.reference}/{wildcards.sample}  
		cd QC_RSeQC/{wildcards.organism}/{wildcards.reference}/{wildcards.sample}  
		{rseqc_bin}/read_quality.py -i ../../../../{input.unsorted_BAM} -o out &> out.read_quality.log
		wait
		{rseqc_bin}/bam_stat.py -i ../../../../{input.unsorted_BAM} &> out.bam_stat.txt & 
		{rseqc_bin}/clipping_profile.py -i ../../../../{input.unsorted_BAM} -s PE -o out &> out.clipping_profile.log & 
		{rseqc_bin}/deletion_profile.py -i ../../../../{input.unsorted_BAM} -l 100 -o out &> out.deletion_profile.log & 
		{rseqc_bin}/FPKM_count.py -i ../../../../{input.sorted_BAM} -o out.only-exonic -r ../../../../{input.RSeQC_reference_bed} -d '1+-,1-+,2++,2--' --only-exonic &> out.only-exonic.FPKM_count.log & 
		{rseqc_bin}/FPKM_count.py -i ../../../../{input.sorted_BAM} -o out.all_reads -r ../../../../{input.RSeQC_reference_bed} -d '1+-,1-+,2++,2--' &> out.all_reads.FPKM_count.log & 
		{rseqc_bin}/geneBody_coverage.py -i ../../../../{input.sorted_BAM} -r ../../../../{input.RSeQC_reference_bed} -o out &> out.geneBody_coverage.log & 
		{rseqc_bin}/infer_experiment.py -i ../../../../{input.unsorted_BAM} -r ../../../../{input.RSeQC_reference_bed} &> out.infer_experiment.txt & 
		{rseqc_bin}/inner_distance.py -i ../../../../{input.unsorted_BAM} -o out -r ../../../../{input.RSeQC_reference_bed} &> out.inner_distance.log & 
		{rseqc_bin}/insertion_profile.py -s PE -i ../../../../{input.unsorted_BAM} -o out &> out.insertion_profile.log &
		{rseqc_bin}/junction_annotation.py -i ../../../../{input.unsorted_BAM} -o out -r ../../../../{input.RSeQC_reference_bed} &> out.junction_annotation.log &
		{rseqc_bin}/junction_saturation.py -i ../../../../{input.unsorted_BAM} -r ../../../../{input.RSeQC_reference_bed} -o out &> out.junction_saturation.log &
		{rseqc_bin}/mismatch_profile.py -l 100 -i ../../../../{input.unsorted_BAM} -o out &> out.mismatch_profile.log &
		{rseqc_bin}/read_distribution.py  -i ../../../../{input.unsorted_BAM} -r ../../../../{input.RSeQC_reference_bed} &> out.read_distribution.txt &
		{rseqc_bin}/read_duplication.py -i ../../../../{input.unsorted_BAM} -o out &> out.read_duplication.log &
		{rseqc_bin}/read_GC.py -i ../../../../{input.unsorted_BAM} -o out &> out.read_GC.log &
		{rseqc_bin}/read_hexamer.py -i ../../../../{input.unsorted_BAM} -r ../../../../{input.RSeQC_reference_bed} -g ../../../../{input.reference_dir}/*.fa &> out.read_hexamer.txt &
		{rseqc_bin}/read_NVC.py -i ../../../../{input.unsorted_BAM} -x -o out &> out.read_NVC.log &
		{rseqc_bin}/RNA_fragment_size.py -r ../../../../{input.RSeQC_reference_bed} -i ../../../../{input.sorted_BAM} &> out.RNA_fragment_size.txt &> out.RNA_fragment_size.log &
		{rseqc_bin}/split_bam.py -i ../../../../{input.unsorted_BAM}  -r ../../../../{input.RSeQC_reference_bed} -o out.split &> out.split.log &
		{rseqc_bin}/tin.py -i ../../../../{input.sorted_BAM} -r ../../../../{input.RSeQC_reference_bed} &> tin.no_subtract_background.txt &
		{rseqc_bin}/tin.py -i ../../../../{input.sorted_BAM} -r ../../../../{input.RSeQC_reference_bed} -s &> tin.subtract_background.txt &
		wait """
