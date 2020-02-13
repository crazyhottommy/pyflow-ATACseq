shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

localrules: all
# localrules will let the rule run locally rather than submitting to cluster
# computing nodes, this is for very small jobs


FILES = json.load(open(config['SAMPLES_JSON']))

CLUSTER = json.load(open(config['CLUSTER_JSON']))

ALL_SAMPLES = sorted(FILES.keys())

ALL_FASTQ   = expand("01seq/{sample}_{read}.fastq.gz", sample = ALL_SAMPLES, read = ["R1", "R2"])
ALL_FASTQC  = expand("02fqc/{sample}_{read}_fastqc.zip", sample = ALL_SAMPLES, read = ["R1", "R2"])
ALL_TRIMMED_FASTQ = expand("03trim/{sample}_{read}.trimmed.fastq.gz", sample = ALL_SAMPLES, read = ["R1", "R2"])

#this is the name sorted bam, not coordinate sorted bam, samtools index only for coordiante sorted bam.
ALL_BAM = expand("04aln/{sample}.sorted.bam", sample = ALL_SAMPLES)

ALL_PHATOM = expand("05phantompeakqual/{sample}_phantom.txt", sample = ALL_SAMPLES)
ALL_FLAGSTAT = expand("04aln/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES)
ALL_CHRM_EXCLUDE_BAM = expand("06aln_exclude_chrM/{sample}_exclude_chrM.sorted.bam", sample = ALL_SAMPLES)
ALL_CHRM_EXCLUDE_BAM_INDEX = expand("06aln_exclude_chrM/{sample}_exclude_chrM.sorted.bam.bai", sample = ALL_SAMPLES)

## downsample from the chrM excluded bam files
ALL_DOWNSAMPLE_BAM = expand("06aln_downsample/{sample}-downsample.sorted.bam", sample = ALL_SAMPLES)
ALL_DOWNSAMPLE_INDEX = expand("06aln_downsample/{sample}-downsample.sorted.bam.bai", sample = ALL_SAMPLES)

ALL_BIGWIG = expand("07bigwig/{sample}.bw", sample = ALL_SAMPLES)
ALL_PEAKS = expand("08peak_macs2/{sample}_macs2_peaks.broadPeak", sample = ALL_SAMPLES)
ALL_NUCLEO = expand("09nucleoATAC/{sample}_nucleoATAC.occpeaks.bed.gz", sample = ALL_SAMPLES)
ALL_QC = ["10multiQC/multiQC_log.html"]
ALL_ATAQV = expand("04aln/{sample}.sorted.bam.ataqv.json", sample = ALL_SAMPLES)
ALL_ATAC_QC = ["11ATAC_qc_html"]


rule all:
	input: ALL_FASTQ + ALL_FASTQC + ALL_TRIMMED_FASTQ + ALL_BAM + ALL_BIGWIG + ALL_FLAGSTAT + ALL_QC + ALL_PEAKS + ALL_PHATOM + \
	ALL_CHRM_EXCLUDE_BAM + ALL_CHRM_EXCLUDE_BAM_INDEX + ALL_NUCLEO + ALL_ATAQV + ALL_ATAC_QC



rule merge_fastqs:
	input:
		r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
		r2 = lambda wildcards: FILES[wildcards.sample]['R2']
	output:
		"01seq/{sample}_R1.fastq.gz", "01seq/{sample}_R2.fastq.gz"
	log: "00log/{sample}_merge_fastq.log"
	params:
		jobname = "{sample}"
	threads: 1
	message: "merging fastqs {input}: {threads} threads"
	shell:
		"""
		gunzip -c {input.r1} | gzip > {output[0]} 2> {log}
		gunzip -c {input.r2} | gzip > {output[1]} 2>> {log}
		"""


rule fastqc:
    input:  "01seq/{sample}_R1.fastq.gz", "01seq/{sample}_R2.fastq.gz"
    output: "02fqc/{sample}_R1_fastqc.zip", "02fqc/{sample}_R2_fastqc.zip"
    log:    "00log/{sample}_fastqc"
    threads: 1
    params : jobname = "{sample}"
    message: "fastqc {input}: {threads}"
    shell:
        """
	# fastqc works fine on .gz file as well
        module load fastqc
        fastqc -o 02fqc -f fastq --noextract {input[0]} {input[1]} 2> {log}
        """

rule trim_adapter:
 	input: "01seq/{sample}_R1.fastq.gz", "01seq/{sample}_R2.fastq.gz"
 	output: "03trim/{sample}_R1.trimmed.fastq.gz", "03trim/{sample}_R2.trimmed.fastq.gz"
 	log: "00log/{sample}_trim_adaptor.log"
 	threads: 1
 	params : jobname = "{sample}"
 	message: "trim_adaptor {input}: {threads}"
 	shell:
 		"""
 		# python2.x
		trim_adapters {input[0]} {input[1]} 2> {log}
		mv 01seq/{wildcards.sample}_R1.trimmed.fastq.gz 03trim/
		mv 01seq/{wildcards.sample}_R2.trimmed.fastq.gz 03trim/
 		"""

## the later step will remove chrM from the bam file and coordinate sort the bam
## so I did not cooridnate sort the bam at this step to save some time.
rule align:
	input: "03trim/{sample}_R1.trimmed.fastq.gz", "03trim/{sample}_R2.trimmed.fastq.gz"
	output: "04aln/{sample}.sorted.bam", "00log/{sample}.align"
	threads: CLUSTER["align"]["cpu"]
	params: jobname = "{sample}"
	message: "aligning {input}: {threads} threads"
	log:
		bowtie2 = "00log/{sample}.align",
		markdup = "00log/{sample}.markdup"
	shell:
		"""
		## samblaster mark duplicates for read id grouped reads. I do not coordinate sort the bam
		bowtie2 --threads 5  -X2000 -x {config[idx_bt2]} -1 {input[0]} -2 {input[1]} 2> {log.bowtie2} \
		| samblaster 2> {log.markdup} \
		| samtools view -Sb - > {output[0]}

		"""

# check number of reads mapped by samtools flagstat
rule flagstat_bam:
    input:  "04aln/{sample}.sorted.bam"
    output: "04aln/{sample}.sorted.bam.flagstat"
    log:    "00log/{sample}.flagstat_bam"
    threads: 1
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """
rule ataqv:
	input: "04aln/{sample}.sorted.bam"
	output: "04aln/{sample}.sorted.bam.ataqv.json"
	log: "00log/{sample}_ataqv.log"
	threads: 1
	params: jobname = "{sample}"
	message: "ataqv quality control for {input}"
	shell:
		"""
		ataqv {config[ataqv_g]} {input} --metrics-file {output} 2> {log}
		"""

rule json_to_html:
	input: expand("04aln/{sample}.sorted.bam.ataqv.json", sample = ALL_SAMPLES)
	output: "11ATAC_qc_html"
	log: "00log/ATAC_qc_html.log"
	threads: 1
	message: "compiling json files to html ATAC-seq QC"
	shell:
		"""
		
		mkarv 11ATAC_qc_html {input}
		"""

## shifting the reads are only critical for TF footprint, for peak calling and making bigwigs, it should be fine using the bams without shifting
# https://sites.google.com/site/atacseqpublic/atac-seq-analysis-methods/offsetmethods
rule remove_chrM_bam:
	input: "04aln/{sample}.sorted.bam"
	output: "06aln_exclude_chrM/{sample}_exclude_chrM.sorted.bam", "06aln_exclude_chrM/{sample}_exclude_chrM.sorted.bam.bai"
	log: "00log/{sample}_exclude_chrM_bam.log"
	threads: 5
	message: "excluding chrM from bam {input} : {threads} threads"
	params: jobname = "{sample}"
	shell:
		"""
		# remove duplicates and reads on chrM, coordinate sort the bam
		# samblaster expects name sorted bamq
		samtools view -h {input} | samblaster --removeDups \
		| grep -v -P '\tchrM\t' \
		| samtools view -Sb -F 4 - \
		| samtools sort -m 2G -@ 5 -T {input}.tmp -o {output[0]}

		samtools index {output[0]}
		"""

rule phantom_peak_qual:
    input: "06aln_exclude_chrM/{sample}_exclude_chrM.sorted.bam"
    output: "05phantompeakqual/{sample}_phantom.txt"
    log: "00log/{sample}_phantompeakqual.log"
    threads: 4
    params: jobname = "{sample}"
    message: "phantompeakqual for {input} : {threads} threads"
    shell:
        """
        Rscript  /scratch/genomic_med/apps/phantompeak/phantompeakqualtools/run_spp_nodups.R -c={input} -savp -rf  -p=4 -odir=05phantompeakqual  -out={output} -tmpdir=05phantompeakqual 2> {log}

        """

## consider how to reuse the rules.
rule flagstat_chrM_exclude_bam:
    input:  "06aln_exclude_chrM/{sample}_exclude_chrM.sorted.bam"
    output: "06aln_exclude_chrM/{sample}_exclude_chrM.sorted.bam.flagstat"
    log:    "00log/{sample}_exclude_chrM_flagstat_bam"
    threads: 1
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

rule down_sample:
    input: "06aln_exclude_chrM/{sample}_exclude_chrM.sorted.bam", "06aln_exclude_chrM/{sample}_exclude_chrM.sorted.bam.bai",
	 	   "06aln_exclude_chrM/{sample}_exclude_chrM.sorted.bam.flagstat"
    output: "06aln_downsample/{sample}-downsample.sorted.bam", "06aln_downsample/{sample}-downsample.sorted.bam.bai"
    log: "00log/{sample}_downsample.log"
    threads: 5
    params: jobname = "{sample}"
    message: "downsampling for {input}"
    run:
        import re
        import subprocess
        with open (input[2], "r") as f:
            # fifth line contains the number of mapped reads
            line = f.readlines()[5]
            match_number = re.match(r'(\d.+) \+.+', line)
			## how many paired reads, roughly total #reads/2
            total_reads = float(match_number.group(1))/2

        target_reads = config["target_reads"] # 15million reads  by default, set up in the config.yaml file
        if total_reads > target_reads:
            down_rate = target_reads/total_reads
        else:
            down_rate = 1

        shell("sambamba view -f bam -t 5 --subsampling-seed=3 -s {rate} {inbam} | samtools sort -m 2G -@ 5 -T {outbam}.tmp > {outbam} 2> {log}".format(rate = down_rate, inbam = input[0], outbam = output[0], log = log))

        shell("samtools index {outbam}".format(outbam = output[0]))


rule make_bigwigs:
    input : "06aln_downsample/{sample}-downsample.sorted.bam", "06aln_downsample/{sample}-downsample.sorted.bam.bai"
    output: "07bigwig/{sample}.bw"
    log: "00log/{sample}.makebw"
    threads: 5
    params: jobname = "{sample}"
    message: "making bigwig for {input} : {threads} threads"
    shell:
        """
    
    	# no window smoothing is done, for paired-end, bamCoverage will extend the length to the fragement length of the paired reads
        bamCoverage -b {input[0]} --ignoreDuplicates --skipNonCoveredRegions --normalizeUsingRPKM -p 5 --extendReads -o {output} 2> {log}

        """
# https://github.com/taoliu/MACS/issues/145
rule call_peaks_macs2:
    input: "06aln_downsample/{sample}-downsample.sorted.bam", "06aln_downsample/{sample}-downsample.sorted.bam.bai"
    output: bed = "08peak_macs2/{sample}_macs2_peaks.broadPeak"
    log: "00log/{sample}_call_peaks_macs2.log"
    params:
        name = "{sample}_macs2",
        jobname = "{sample}"
    message: "call_peaks macs2 {input}: {threads} threads"
    shell:
        """
       
       ## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
        macs2 callpeak -t {input[0]} \
            --keep-dup all -f BAMPE -g {config[macs2_g]} \
            --outdir 08peak_macs2 -n {params.name} -p {config[macs2_pvalue]} \
            --broad --broad-cutoff {config[macs2_pvalue_broad]} &> {log}
        """

rule multiQC:
    input :
        expand("00log/{sample}.align", sample = ALL_SAMPLES),
        expand("04aln/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES),
        expand("02fqc/{sample}_{read}_fastqc.zip", sample = ALL_SAMPLES, read = ["R1", "R2"])
    output: "10multiQC/multiQC_log.html"
    log: "00log/multiqc.log"
    message: "multiqc for all logs"
    shell:
        """
        multiqc 02fqc 04aln 00log -o 10multiQC -d -f -v -n multiQC_log 2> {log}

        """

## extend the broad peak a bit for nucelosome analysis by nuceloATAC
rule make_bed_nucleoATAC:
	input: "08peak_macs2/{sample}_macs2_peaks.broadPeak"
	output: "09nucleoATAC/{sample}_nucleo.bed"
	log: "00log/{sample}_make_nucleoATAC_bed.log"
	threads: 1
	message: "making bed for nucleoATAC from {input}"
	params: jobname= "{sample}"
	shell:
		"""
		cat {input} | bedtools slop -b 200 -g {config[genome_size]} | sort -k1,1 -k2,2n | bedtools merge > {output} 2> {log}
		"""

## nucleoATAC works on the non-shifted bam, and shift the reads internally!
# https://github.com/GreenleafLab/NucleoATAC/issues/58
rule nucleoATAC:
	input: "06aln_downsample/{sample}-downsample.sorted.bam", "06aln_downsample/{sample}-downsample.sorted.bam.bai", "09nucleoATAC/{sample}_nucleo.bed"
	output: "09nucleoATAC/{sample}_nucleoATAC.occpeaks.bed.gz"
	log: "00log/{sample}_nucleoATAC.log"
	threads: 5
	message: "calling nucleosome by nucleoATAC for {input} : {threads} threads"
	params:
		jobname = "{sample}",
		outputdir = os.path.dirname(srcdir("00log"))
	shell:
		"""
		
		cd 09nucleoATAC
		nucleoatac run --bed {params.outputdir}/{input[2]} --bam {params.outputdir}/{input[0]} --cores 5 --fasta {config[genome_fasta]} --out {wildcards.sample} 2> {params.outputdir}/{log}
		"""
