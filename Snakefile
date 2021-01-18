import os
import numpy as np
from csv import DictReader

#configfile: "config.yaml"
path_hilo = os.getcwd() + "/" # get absolute path to this hilo git directory on the local machine

# wildcards
wildcard_constraints:
    ID = "[A-Za-z0-9]+",
    POP = "pop[0-9]+|allopatric_maize|allopatric_maize_subsample[0-9]+|sympatric_maize|sympatric_mexicana|allopatric_mexicana", # for f4 stats allopatric_maize needs to be treated like a pop
    BOOT = "[0-9]+", # bootstrap
    YESNO = "yes|no",
    Ne = "[0-9]+",
    POSNEG = "pos|neg",
    SIG = "fdr05|perc02|p05", # outlier significance cutoffs: 5% FDR, 2% empirical cutoff, p = 5% (5% percentile from null simulations)
    STAT = "meanAnc|lmElev", # statistics defining outliers
    GROUP = "sympatric_maize|sympatric_mexicana|allopatric_maize|allopatric_mexicana|parv",
    ZEA = "maize|mexicana",
    ALLO_MEX = "allopatric_maize|pop22", # used in f4s. pop22 is Amecameca (at 2467m, the highest elevation of 3 allopatric pops)
    SUBSAMPLE = "[0-9]+",
    r = "[0-9]+", # numeric only
    FEATURE = "r|cd|frac", # recombination rate cM/Mb (r), coding bp/cM (cd), or frac coding bp (frac)
    COVERAGE = "ALL|Over0.5x" # we only estimate local ancestry for individuals with >0.5x mean coverage. whereas allopatric pops and global ancestry estimates we used a less stringent cutoff (individuals with >0.1x coverage)

# reference genome and associated files
#ref = "/home/ecalfee/hilo/data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
ref = path_hilo + "data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
fai = ref + ".fai"
ref_chr = path_hilo + "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
ref_anno = path_hilo + "data/refMaize/geneAnnotations/Zea_mays.B73_RefGen_v4.41.chr.gff3"
trip_anc = "filtered_bams/results/SRR7758238/TRIP.fa.gz"

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


# all bams included in combined sample (est_coverage > 0.05x)
prefix_all = "HILO_MAIZE55"
#prefix_all = "TEST"

# list of all included bams and ids (over minimum 0.05x coverage)
with open("samples/" + prefix_all + "_bams.list") as f:
    all_bams = f.read().splitlines()
with open("samples/" + prefix_all + "_ids.list") as f:
    all_ids = f.read().splitlines()
# add bams for outgroup tripsacum
with open("samples/ALL_byPop/trip_bams.list") as f:
    trip_bams = f.read().splitlines()
# add bams for parviglumis
with open("samples/ALL_byPop/parv_bams.list") as f:
    parv_bams = f.read().splitlines()
all_bams_parv_trip = all_bams + parv_bams + trip_bams

# samples with local ancesty calls (sympatric and over 0.5x coverage)
#with open("samples/Over0.5x_byPop/sympatric_maize_ids.list") as f:
#        symp_maize_ids = f.read().splitlines()
#with open("samples/Over0.5x_byPop/sympatric_mexicana_ids.list") as f:
#        symp_mex_ids = f.read().splitlines()
#symp_ids = symp_maize_ids + symp_mex_ids

# groups
groups = ["sympatric_maize", "sympatric_mexicana", "allopatric_maize", "allopatric_mexicana"]
allo_groups = ["allopatric_maize", "allopatric_mexicana"]
zea = ["maize", "mexicana"]
Nes = [1000, 10000, 100000]

# sympatric populations (maize and mexicana)
symp_pops = ["pop18", "pop19", "pop21", "pop23", "pop24", "pop25", "pop26",
"pop27", "pop28", "pop29", "pop30", "pop31", "pop34", "pop35",
"pop360", "pop361", "pop362", "pop363", "pop365", "pop366", "pop367",
"pop368", "pop369", "pop370", "pop371", "pop372", "pop373", "pop374"]
# allopatric mexicana pops
allo_mex_pops = ["pop20", "pop22", "pop33"]
# includes all pops sequenced by this study (allopatric maize, parviglumis & tripsacum published by other studies)
hilo_pops = allo_mex_pops + symp_pops

# create a dictionary that has one entry for each sympatric population
# containing a list of included id's for samples over 0.5x coverage (e.g. local ancestry inference)
# for each sympatric population.
symp_dict = {}
for pop in symp_pops:
    with open("samples/Over0.5x_byPop/" + pop + "_ids.list") as f:
        symp_dict[pop] = f.read().splitlines()


# make a dictionary of 5Mb regions across the genome
regions_dict = {}
with open("data/refMaize/divide_5Mb/ALL_regions.list") as f:
    for line in f:
        row = line.split("\t")
        # one key, value dictionary entry per region, e.g. regions_dict[6] = 1:30000000-34999999 for region 6
        regions_dict["region_" + row[3]] = row[0] + ":" + row[1] + "-" + row[2]


# snakemake sub-workflows
# note: commenting out some workflows that are already completed makes DAG a lot faster!
#include: "filtered_bams/Snakefile"
#include: "variant_sites/Snakefile"
#include: "global_ancestry/Snakefile"
include: "local_ancestry/Snakefile"
include: "ancestry_by_r/Snakefile"
#include: "ZAnc/Snakefile"
include: "diversity/Snakefile"

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
        # block bootstrap of ancestry (NGSAdmix) ~ recombination rate quintile
        #expand("ancestry_by_r/results/BED_1cM/{WINDOW}.bed", WINDOW = windows_1cM),
        #expand("ancestry_by_r/esults/bootstrap_1cM/" + prefix_all + "/r5_{r}/boot{BOOT}.list",
        #r = [1, 2, 3, 4, 5], BOOT = list(range(0,101))),
        #expand("ancestry_by_r/results/GL_1cM/" + prefix_all + "/{WINDOW}.beagle.gz", WINDOW = windows_1cM),
        #expand("ancestry_by_r/results/bootstrap_1cM/" + prefix_all + "/r5_{r}/boot{BOOT}.beagle.gz",
        #r = [1, 2, 3, 4, 5], BOOT = list(range(0,101))),
        #expand("ancestry_by_r/results/bootstrap_1cM/" + prefix_all + "/r5_{r}/K2/boot{BOOT}.anc",
        #r = [1, 2, 3, 4, 5], BOOT = list(range(0,101))),
        #expand("ancestry_by_r/results/bootstrap_1cM/" + prefix_all + "/{FEATURE}5_K2.Rdata", FEATURE = ["r", "cd"]),
        "ancestry_by_r/plots/K2_by_r_bootstrap_sympatric_only.png",
        "ancestry_by_r/plots/K2_by_r_bootstrap_sympatric_and_allopatric.png",
        "ancestry_by_r/plots/K2_by_cd_bootstrap_sympatric_only.png",
        "ancestry_by_r/plots/K2_by_cd_bootstrap_sympatric_and_allopatric.png",
        "ancestry_by_r/plots/K2_by_r_bootstrap_lm_elevation_facet_r.png",
        "ancestry_by_r/plots/K2_by_r_bootstrap_lm_elevation_color_elev.png",
        "ancestry_by_r/tables/elev_r_interaction.tex",
        "ancestry_by_r/plots/local_anc_by_r_sympatric_only.png",
        "ancestry_by_r/plots/local_anc_by_cd_sympatric_only.png",
        "ancestry_by_r/plots/boot_cor_local_anc_bp.png",
        "ancestry_by_r/plots/boot_cor_local_anc_r5.png",
        "ancestry_by_r/tables/pearsons_rho_ngsadmix.tex",
        expand("ancestry_by_r/tables/tbl_pearsons_rho_f4_{POP}_{ALLO_MEX}.tex", POP = ["sympatric_maize", "sympatric_mexicana"], ALLO_MEX = "pop22"),
        # local ancestry inference
        #expand("local_ancestry/results/alloFreqs/" + prefix_all + "/{GROUP}/{REGION}.mafs.gz", GROUP=allo_groups, REGION=list(regions_dict.keys())),
        #"local_ancestry/results/thinnedSNPs/" + prefix_all + "/whole_genome.var.sites",
        "local_ancestry/results/thinnedSNPs/" + prefix_all + "/whole_genome.bed",
        #expand("local_ancestry/results/countsMajMin/" + prefix_all + "/{ID}.counts.txt", ID = all_ids),
        #expand("local_ancestry/results/ancestry_hmm/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/anc/{POP}.anc.freq",
        expand("local_ancestry/results/ancestry_hmm/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/MAP/{POP}.anc.ind",
        Ne = 10000, YESNO = "yes", POP = symp_pops),
        expand("local_ancestry/results/ancestry_hmm/" + prefix_all + "/Ne{Ne}_yesBoot/anc/{ZEA}.combined.anc.bed", ZEA = zea, Ne = Nes),
        #expand("ancestry_by_r/results/local_anc_1cM/" + prefix_all + "/Ne10000_yesBoot/{POP}.bed", POP = symp_pops),
        expand("ancestry_by_r/results/local_anc_1cM/" + prefix_all + "/Ne{Ne}_yesBoot/{POP}.anc.wind", POP = symp_pops, Ne = Nes),
        expand("local_ancestry/results/ancestry_hmm/" + prefix_all + "/Ne{Ne}_yesBoot/{POP}.times", POP = symp_pops, Ne = Nes),
        expand("local_ancestry/plots/admix_times_Ne{Ne}_{YESNO}Boot.png", Ne = Nes, YESNO = "yes"),
        expand("ZAnc/results/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/{ZEA}.MVN.RData", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("ZAnc/results/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/{ZEA}.lmElev.fit.RData", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("ZAnc/plots/Ne{Ne}_{YESNO}Boot/{ZEA}_hist_outlier_peaks.png", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("ZAnc/plots/Ne{Ne}_{YESNO}Boot/{ZEA}_slope_elev.png", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("ZAnc/plots/Ne{Ne}_{YESNO}Boot/{ZEA}_mean_anc.png", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("ZAnc/results/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/{ZEA}.zAnc.fdr.RData", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("ZAnc/results/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/{ZEA}.zAnc.fit.RData", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("local_ancestry/results/ancestry_hmm/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/HOMOZYG/{ZEA}/bams/{POP}.completed", POP = symp_pops, Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("local_ancestry/results/ancestry_hmm/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/HOMOZYG/{ZEA}/bams/{POP}_bams.list", POP = symp_pops, Ne = 10000, ZEA = zea, YESNO = "yes"),
        trip_anc,
        expand("ancestry_by_r/results/f4/{POP}.Dstats.Observed.txt", ALLO_MEX = "pop22", POP = ["sympatric_maize", "sympatric_mexicana", "pop22"]),
        expand("ancestry_by_r/results/f4/{POP}.f4", POP = ["sympatric_maize", "sympatric_mexicana", "pop22"]),
        expand("ancestry_by_r/plots/f4_{POP}_{ALLO_MEX}_byr5.png", ALLO_MEX = "pop22", POP = ["sympatric_maize", "sympatric_mexicana"]),
        expand("ZAnc/results/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/flowering_time_genes_v4.plus10kb.{ZEA}_{POSNEG}_{STAT}_outliers.{SIG}.{SUFFIX}", ZEA = zea, Ne = 10000, YESNO = "yes", POSNEG = ["pos", "neg"], STAT = ["meanAnc", "lmElev"], SIG = ["fdr05", "perc02", "p05"], SUFFIX = ["bed", "counts"]),
        expand("diversity/results/pi/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/HOMOZYG/{ZEA}/{POP}.thetas.gz", Ne = 100000, YESNO = "yes", ZEA = zea, POP = symp_pops)
    params:
        p = "med2"
    resources:
        time_min = 30,
        mem = 2

## some: alternative to all for running part of the pipeline (e.g. testing or pipeline incomplete)
rule some:
    input:

    params:
        p = "med2"
    resources:
        time_min = 15,
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
