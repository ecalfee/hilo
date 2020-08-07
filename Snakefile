import os
import numpy as np
from csv import DictReader

#configfile: "config.yaml"
path_hilo = os.getcwd() + "/" # get absolute path to this hilo git directory on the local machine

# reference genome
#ref = "/home/ecalfee/hilo/data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
ref = path_hilo + "data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
fai = ref + ".fai"
ref_chr = path_hilo + "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"

# recombination map
rmap = path_hilo + "data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt"
rmap_ext = path_hilo + "data/linkage_map/ogut_fifthcM_map_agpv4_EXTENDED.txt"

# 1cM genomic windows for whole genome (1520 windows)
windows_1cM = []
#with open("ancestry_by_r/results/map_pos_1cM_windows.txt", "r") as read_obj:
#    csv_dict_reader = DictReader(read_obj, delimiter = "\t")
#    for row in csv_dict_reader:
#        windows_1cM.append(row["window"])
for i in range(1, 1521):
    windows_1cM.append("W" + str(i))

# all fastq/samples sequenced
#prefix_bams = "April2020"
#prefix_bams = "TEST2"
prefix_bams = "Combined"

# all bams included in combined sample (est_coverage > 0.05x)
prefix_all = "HILO_MAIZE55"
#prefix_all = "TEST"

# groups
groups = ["sympatric_maize", "sympatric_mexicana", "allopatric_maize", "allopatric_mexicana"]
allo_groups = ["allopatric_maize", "allopatric_mexicana"]

# sympatric populations
sympatric_pops = ["pop18", "pop19", "pop21", "pop23", "pop24", "pop25", "pop26",
"pop27", "pop28", "pop29", "pop30", "pop31", "pop34", "pop35",
"pop360", "pop361", "pop362", "pop363", "pop365", "pop366", "pop367",
"pop368", "pop369", "pop370", "pop371", "pop372", "pop373", "pop374"]

# make a dictionary of 5Mb regions across the genome
regions_dict = {}
with open("data/refMaize/divide_5Mb/ALL_regions.list") as f:
    for line in f:
        row = line.split("\t")
        # one key, value dictionary entry per region, e.g. regions_dict[6] = 1:30000000-34999999 for region 6
        regions_dict["region_" + row[3]] = row[0] + ":" + row[1] + "-" + row[2]


# snakemake sub-workflows
#include: "filtered_bams/Snakefile"
#include: "variant_sites/Snakefile"
include: "global_ancestry/Snakefile"
include: "local_ancestry/Snakefile"
#include: "ancestry_by_r/Snakefile"

## all:  main rule to run all workflows
rule all:
    input:
        # SNP set
        expand("variant_sites/results/" + prefix_all + "/{REGION}.rpos",
                REGION=list(regions_dict.keys())),
        # bam metrics files
        "filtered_bams/metrics/fastQC/multiqc/multiqc_report.html",
        "filtered_bams/metrics/fastQC_trimmed/multiqc/multiqc_report.html",
        "filtered_bams/metrics/picard/multiqc/multiqc_report.html",
        "filtered_bams/metrics/flagstat/multiqc/multiqc_report.html",
        # all bams
        #expand(["filtered_bams/merged_bams/{ID}.sort.dedup.bam",
        #        "filtered_bams/merged_bams/{ID}.sort.dedup.bam.bai"],
        #        zip, ID=list(merge_dict.keys())),
        # global ancestry analysis: thinned GL, PCA and NGSAdmix
        "global_ancestry/plots/pca.png",
        "global_ancestry/plots/lm_mexicana_by_pop_elevation_K2.png",
        #"global_ancestry/results/thinnedSNPs/" + prefix_all + "/whole_genome.beagle.gz",
        #"global_ancestry/results/PCA/" + prefix_all + "/whole_genome.cov",
        #"global_ancestry/results/NGSAdmix/" + prefix_all + "/K2.qopt",
        "filtered_bams/plots/p_seq_counts.png",
        expand("samples/{SUBSET}_byPop/{GROUP}_{LIST_TYPE}.list", GROUP=groups, LIST_TYPE=["ids", "bams"], SUBSET=["ALL", "Over0.5x"]),
        expand("samples/{SUBSET}_byPop/{POP}_{LIST_TYPE}.list", POP=sympatric_pops, LIST_TYPE=["ids", "bams"], SUBSET=["ALL", "Over0.5x"]),
        expand("variant_sites/results/popFreq/{GROUP}/{REGION}.mafs.gz", GROUP=allo_groups, REGION=list(regions_dict.keys())),
        # block bootstrap of ancestry (NGSAdmix) ~ recombination rate quintile
        #expand("ancestry_by_r/results/BED_1cM/{WINDOW}.bed", WINDOW = windows_1cM),
        #expand("ancestry_by_r/esults/bootstrap_1cM/" + prefix_all + "/r5_recomb{r}/boot{BOOT}.list",
        #r = [1, 2, 3, 4, 5], BOOT = list(range(0,101))),
        #expand("ancestry_by_r/results/GL_1cM/" + prefix_all + "/{WINDOW}.beagle.gz", WINDOW = windows_1cM),
        #expand("ancestry_by_r/results/bootstrap_1cM/" + prefix_all + "/r5_recomb{r}/boot{BOOT}.beagle.gz",
        #r = [1, 2, 3, 4, 5], BOOT = list(range(0,101))),
        #expand("ancestry_by_r/results/bootstrap_1cM/" + prefix_all + "/r5_recomb{r}/K2/boot{BOOT}.anc",
        #r = [1, 2, 3, 4, 5], BOOT = list(range(0,101))),
        "ancestry_by_r/plots/K2_by_r_bootstrap_sympatric_only.png",
        "ancestry_by_r/plots/K2_by_r_bootstrap_sympatric_and_allopatric.png",
        "ancestry_by_r/plots/K2_by_r_bootstrap_lm_elevation_facet_r.png",
        "ancestry_by_r/plots/K2_by_r_bootstrap_lm_elevation_color_elev.png"
    params:
        p = "med2"
    resources:
        time_min = 15,
        mem = 2



## some: alternative to all for running part of the pipeline (e.g. testing or pipeline incomplete)
rule some:
    input:
        # SNP set for 1 scaffold
        #"variant_sites/results/" + prefix_all + "/region_1.rpos",
        #"variant_sites/results/" + prefix_all + "/region_1.var.sites"
        "global_ancestry/results/thinnedSNPs/" + prefix_all + "/whole_genome.beagle.gz",
        expand("variant_sites/results/popFreq/{GROUP}/{REGION}.mafs.gz", GROUP=["allopatric_mexicana"], REGION=["region_120"])
    params:
        p = "med2"
    resources:
        time_min = 5,
        mem = 2

## test: for running test files
rule test:
    input:
        "test/whole_genome.beagle.gz",
        #"test/whole_genome.cov"
        "test/K2.qopt"
    params:
        p = "med2"
    resources:
        time_min = 5,
        mem = 2
