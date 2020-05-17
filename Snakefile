rule all:
    input:
        "results/TEST/a.sort.dedup.baq.bam.bai",
        "results/TEST/b.sort.dedup.baq.bam.bai"
        #"plots/quals.svg"

ref="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"

# align fastq to maize reference AGPv4 official release using bwa mem
# note: reference genome needs to already be indexed, e.g. bwa index ref.fa
rule index_ref:
    input:
        "{ref}"
    output:
        "{ref}.fai"
    #log:
    #    "logs/index_ref_bwa.log"
    #conda:
    #    "environment.yaml"
    #envmodules:
    #    "bio/bwa/0.7.9",
    #    "bio/samtools/1.9"
    shell:
        "bwa index {input}"


# LIBRARY WILDCARD?
rule bwa_map:
    input:
        "{ref}",
        "{ref}.fai",
        "../data/HILO_raw_reads/{LANE}/{ID}_1.fq.gz", # fastq1
        "../data/HILO_raw_reads/{LANE}/{ID}_2.fq.gz" # fastq2
    output:
        #temp("results/{LANE}/{ID}.bam") # unsorted tmp bams are deleted once consuming jobs are completed
        "results/{LANE}/{ID}.bam"
    params:
        rg="@RG\tID:{LANE}\tSM:{ID}\tPL:ILLUMINA\tLB:{LIBRARY}\tPU:{ID}.{LANE}"
    #threads: 16
    #log:
    #    "logs/bwa_map/{ID}_{LANE}.log"
    #benchmark:# record memory usage and time
    #   "benchmarks/{ID}_{LANE}.bwa_map.txt"
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} -v 3 {input} | "
        "samtools view -Shu -) > {output} 2> {log}"


#4 threads. tmp directory
rule samtools_sort:
    input:
        "results/{LANE}/{ID}.bam"
    output:
        "results/{LANE}/{ID}.sort.bam"
    #log:
    #    "logs/samtools_sort/{ID}_{LANE}.log"
    shell:
        "mkdir -p /scratch/ecalfee/sorting_reads/{wildcards.ID}_{wildcards.LANE};"
        "samtools sort -m {mem} -@ {threads} "
        "-T /scratch/ecalfee/sorting_reads/{wildcards.ID}_{wildcards.LANE} "
        "-O bam {input} > {output}"

rule remove_duplicates:
    input:
        "{ref}",
        "{ref}.fai",
        "results/{LANE}/{ID}.sort.bam"
    output:
        bam="results/{LANE}/{ID}.dedup.baq.bam",
        metrics="metrics/{LANE}/{ID}.metrics.txt"
    #log:
    #    "logs/remove_duplicates/{ID}_{LANE}.log"
    shell:
        "mkdir -p /scratch/ecalfee/dedup_reads/{wildcards.ID}_{wildcards.LANE};"
        "java -Xmx6g -jar ${PICARD}/picard.jar MarkDuplicates "
        "INPUT= {input} OUTPUT=/dev/stdout QUIET=true "
        "TMP_DIR=/scratch/ecalfee/dedup_reads/{wildcards.ID}_{wildcards.LANE} "
        "METRICS_FILE={results.metrics} | "
        "samtools calmd -SArE --reference {REF} - | "
        "samtools view -bS -q 30 - > {results.bam}"

rule samtools_index:
    input:
        "results/{LANE}/{ID}.sort.dedup.baq.bam"
    output:
        "results/{LANE}/{ID}.sort.dedup.baq.bam.bai"
    shell:
        "samtools index {input}"

# now merge any bams for IDs with multiple sequencing runs
# TO-DO:
# load hilo ID's from a file into an array
# load lanes into an array also
# make 2 test fastq sets _1 and _2 (very small). put in data/HILO_raw_reads/
# put this snakemake file within filtered_bams and run from that working directory.
# make snake_logs directory within filtered_bams
# activate conda environment. try to run the snakemake. make dag.
# figure out if I can activate environment before calling slurm (unlikely) or test activating with --use-conda for each rule or universally, e.g. in submit.json
