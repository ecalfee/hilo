import os
import numpy as np

#configfile: "config.yaml"
path_hilo = os.getcwd() + "/" # get absolute path to this hilo git directory on the local machine

#ref = "/home/ecalfee/hilo/data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
ref = path_hilo + "data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
fai = ref + ".fai"
ref_chr = path_hilo + "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"

# all fastq/samples sequenced
#prefix_bams = "April2020"
#prefix_bams = "TEST2"
prefix_bams = "Combined"

# all bams included in combined sample (est_coverage > 0.05x)
prefix_all = "HILO_MAIZE55"
#prefix_all = "TEST"

# all bams for local ancestry inference (est_coverage > 0.5x)

# snakemake sub-workflows
include: "filtered_bams/Snakefile"
include: "variant_sites/Snakefile"



rule all:
    input:
        # SNP set
        expand(["variant_sites/results/" + prefix_all + "/region_{REGION_N}.mafs.gz",
                "variant_sites/results/" + prefix_all + "/region_{REGION_N}.beagle.gz"], # multiext("some/plot", ".pdf", ".svg", ".png")
                REGION_N=list(regions_dict.keys())),
        # bam metrics files
        "filtered_bams/metrics/fastQC/multiqc/multiqc_report.html",
        "filtered_bams/metrics/fastQC_trimmed/multiqc/multiqc_report.html",
        "filtered_bams/metrics/picard/multiqc/multiqc_report.html",
        "filtered_bams/metrics/flagstat/multiqc/multiqc_report.html",
        # all bams
        expand(["filtered_bams/merged_bams/{ID}.sort.dedup.bam",
                "filtered_bams/merged_bams/{ID}.sort.dedup.bam.bai"],
                zip, ID=list(merge_dict.keys()))
    params:
        p = "med2"
    resources:
        time_min = 5,
        mem = 2



# alternative to all for running part of the pipeline (e.g. testing or pipeline incomplete)
rule some:
    input:
        # SNP set for 1 scaffold
        "variant_sites/results/" + prefix_all + "/region_1.mafs.gz"
    params:
        p = "med2"
    resources:
        time_min = 5,
        mem = 2
