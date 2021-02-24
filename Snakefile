import os
import numpy as np
from csv import DictReader
import itertools

#configfile: "config.yaml"
path_hilo = os.getcwd() + "/" # get absolute path to this hilo git directory on the local machine

# wildcards
wildcard_constraints:
    ID = "[A-Za-z0-9]+",
    POP = "pop[0-9]+|allopatric_maize|allopatric_maize_subsample[0-9]+|sympatric_maize|sympatric_mexicana|allopatric_mexicana", # for f4 stats allopatric_maize needs to be treated like a pop
    NOTPOP = "not[0-9]+", # eg. not360 is all maize EXCEPT population 360
    POP1 = "pop[0-9]+|not[0-9]+", # pop 1 and pop 2 for pairwise Fst calculations and saf/sfs. Can be 1 population (e.g. pop360) or all individuals in that subspecies except 1 population (e.g. not360)
    POP2 = "pop[0-9]+|not[0-9]+",
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
    COVERAGE = "ALL|Over0.5x", # we only estimate local ancestry for individuals with >0.5x mean coverage. whereas allopatric pops and global ancestry estimates we used a less stringent cutoff (individuals with >0.1x coverage)
    WIN = "[0-9]+", # window size
    STEP = "[0-9]+", # window step size (non-overlapping windows if STEP = WIN)
    n = "[0-9]+"

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
# (in order of elevation)
symp_mexicana_pops = ["pop29", "pop27", "pop19", "pop28", "pop18", "pop21", "pop34", "pop26", "pop35", "pop25", "pop23", "pop30", "pop24", "pop31"]
symp_maize_pops = ["pop371", "pop369", "pop361", "pop370", "pop360", "pop363", "pop362", "pop368", "pop374", "pop367", "pop365", "pop372", "pop366", "pop373"]
symp_pops = symp_mexicana_pops + symp_maize_pops
# allopatric mexicana pops
allo_mex_pops = ["pop20", "pop22", "pop33"]
# includes all pops sequenced by this study (allopatric maize, parviglumis & tripsacum published by other studies)
hilo_pops = allo_mex_pops + symp_pops

symp_mexicana_nots = ["not29", "not27", "not19", "not28", "not18", "not21", "not34", "not26", "not35", "not25", "not23", "not30", "not24", "not31"]
symp_maize_nots = ["not371", "not369", "not361", "not370", "not360", "not363", "not362", "not368", "not374", "not367", "not365", "not372", "not366", "not373"]
symp_nots = symp_mexicana_nots + symp_maize_nots

# get all unique pairs of populations (all combinations with no same-same pairs):
symp_maize_pairs = []
for i in itertools.combinations(symp_maize_pops, 2):
    symp_maize_pairs.append(i[0] + "." + i[1])
symp_mexicana_pairs = []
for i in itertools.combinations(symp_mexicana_pops, 2):
    symp_mexicana_pairs.append(i[0] + "." + i[1])


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
#include: "local_ancestry/Snakefile"
#include: "ancestry_by_r/Snakefile"
include: "ZAnc/Snakefile"
include: "diversity/Snakefile"
#include: "map/Snakefile"

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
        "ancestry_by_r/tables/elev_r_interaction_5.tex",
        "ancestry_by_r/tables/spearmans_rho_ngsadmix.tex",
        "ancestry_by_r/plots/local_anc_by_r_quintiles.png",
        "ancestry_by_r/plots/local_anc_by_cd_quintiles.png",
        "ancestry_by_r/plots/local_anc_by_r_continuous.png",
        "ancestry_by_r/tables/spearmans_rho_local_ancestry.tex",
        expand("ancestry_by_r/tables/spearmans_rho_f4_{POP}_{ALLO_MEX}.tex", POP = ["sympatric_maize", "sympatric_mexicana"], ALLO_MEX = "pop22"),
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
        expand("local_ancestry/results/admix_times_Ne{Ne}_{YESNO}Boot.{SUFFIX}", Ne = Nes, YESNO = "yes", SUFFIX = ["txt", "RDS"]),
        expand("local_ancestry/plots/admix_times_Ne{Ne}_{YESNO}Boot.png", Ne = Nes, YESNO = "yes"),
        expand("ZAnc/results/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/{ZEA}.MVN.RData", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("ZAnc/results/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/{ZEA}.lmElev.fit.RData", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("ZAnc/plots/Ne{Ne}_{YESNO}Boot/mex_maize_hist_outlier_peaks.png", Ne = 10000, YESNO = "yes"),
        expand("ZAnc/plots/Ne{Ne}_{YESNO}Boot/{ZEA}_slope_elev.png", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("ZAnc/plots/Ne{Ne}_{YESNO}Boot/{ZEA}_mean_anc.png", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("ZAnc/plots/Ne{Ne}_{YESNO}Boot/multi_maize_mexicana_genome_scan.png", Ne = 10000, YESNO = "yes"),
        expand("ZAnc/plots/Ne{Ne}_{YESNO}Boot/network_peak_sharing_data_only.png", Ne = 10000, YESNO = "yes"),
        expand("ZAnc/results/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/{ZEA}.zAnc.fdr.RData", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("ZAnc/results/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/{ZEA}.zAnc.fit.RData", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("ZAnc/plots/Ne{Ne}_{YESNO}Boot/QQ.png", Ne = 10000, YESNO = "yes"),
        expand("ZAnc/plots/Ne{Ne}_{YESNO}Boot/combmatrix_peak_sharing_{ZEA}.png", Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("local_ancestry/results/ancestry_hmm/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/HOMOZYG/{ZEA}/bams/{POP}.completed", POP = symp_pops, Ne = 10000, ZEA = zea, YESNO = "yes"),
        expand("local_ancestry/results/ancestry_hmm/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/HOMOZYG/{ZEA}/bams/{POP}_bams.list", POP = symp_pops, Ne = 10000, ZEA = zea, YESNO = "yes"),
        trip_anc,
        expand("ancestry_by_r/results/f4/{POP}.Dstats.Observed.txt", ALLO_MEX = "pop22", POP = ["sympatric_maize", "sympatric_mexicana", "pop22"]),
        expand("ancestry_by_r/results/f4/{POP}.f4", POP = ["sympatric_maize", "sympatric_mexicana", "pop22"]),
        expand("ancestry_by_r/plots/f4_{POP}_{ALLO_MEX}_byr5.png", ALLO_MEX = "pop22", POP = ["sympatric_maize", "sympatric_mexicana"]),
        expand("ZAnc/results/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/flowering_time_genes_v4.plus20kb.{ZEA}_{POSNEG}_{STAT}_outliers.{SIG}.{SUFFIX}", ZEA = zea, Ne = 10000, YESNO = "yes", POSNEG = ["pos", "neg"], STAT = ["meanAnc", "lmElev"], SIG = ["fdr05", "perc02", "p05"], SUFFIX = ["bed", "counts"]),
        expand("ZAnc/results/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/flowering_time_genes_v4.plus20kb.overlap.summary_overlap_outliers.txt", Ne = 10000, YESNO = "yes"),
        expand("diversity/results/pi/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/HOMOZYG/{ZEA}/{POP}.thetas.gz", Ne = 10000, YESNO = "yes", ZEA = zea, POP = symp_pops),
        expand("diversity/results/pi/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/HOMOZYG/{ZEA}/{POP}.pi.windows.{WIN}.{STEP}.pestPG", WIN = 5000, STEP = 5000, Ne = 10000, YESNO = "yes", ZEA = zea, POP = symp_pops),
        expand("diversity/results/pi/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/HOMOZYG/{ZEA}/{POP}.pi.allChr.pestPG", Ne = 10000, YESNO = "yes", ZEA = zea, POP = symp_pops),
        expand("diversity/results/pi/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/HOMOZYG/{ZEA}/{POP}.pi.outliers_vs_not.{WIN}.{STEP}.txt", WIN = 5000, STEP = 5000, Ne = 10000, YESNO = "yes", ZEA = zea, POP = symp_pops),
        "map/plots/mexico_lines_elev_teo_color.png",
        "map/plots/mexico_lines_elev_teo_black.png",
        "ZAnc/tables/" + prefix_all + "/Ne10000_yesBoot/genes_mapped_to_outliers.tex",
        "global_ancestry/plots/global_anc_multi.png"
    params:
        p = "med2"
    resources:
        time_min = 30,
        mem = 2

## fst: alternative to all for running part of the pipeline (e.g. testing or pipeline incomplete)
rule fst:
    input:
        # get fst between sympatric maize and mexicana pairs for both high confidence maize and mexicana ancestry.
        expand("diversity/results/fst/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/mexicana/{POP1}.{POP2}.fst.allChr.txt", POP1 = symp_maize_pops, POP2 = symp_mexicana_pops), # fst between subspecies
        expand("diversity/results/fst/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/maize/{POP1}.{POP2}.fst.allChr.txt", POP1 = symp_mexicana_pops, POP2 = symp_maize_pops),
        expand("diversity/results/fst/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/mexicana/{POP_PAIR}.fst.allChr.txt", POP_PAIR = symp_maize_pairs), # within subspecies, fst for introgressed ancestry
        expand("diversity/results/fst/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/maize/{POP_PAIR}.fst.allChr.txt", POP_PAIR = symp_mexicana_pairs),
        expand("diversity/results/fst/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/mexicana/{POP_PAIR}.fst.allChr.txt", POP_PAIR = symp_mexicana_pairs), # within subspecies, fst for native ancestry
        expand("diversity/results/fst/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/maize/{POP_PAIR}.fst.allChr.txt", POP_PAIR = symp_maize_pairs),
        "diversity/results/fst/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/summary_pop_pairs_fst.allChr.txt",
        expand("diversity/plots/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/fst_within_maize_or_mexicana_ancestry_genomewide_heatmap_both.png", Ne = 10000, YESNO = "yes"),
    params:
        p = "med2"
    resources:
        time_min = 30,
        mem = 2

rule fst_mexicana_anc_outliers:
    input:
        expand("diversity/results/pi/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/mexicana/{POP}.{n}pop.outliers{POP}.pi.allChr.pestPG", POP = symp_maize_pops, n = [1, 4]), # introgressed mexicana tracts within sympatric maize
        expand("diversity/results/pi/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/mexicana/{POP2}.1pop.outliers{POP1}.pi.allChr.pestPG", zip, POP1 = symp_maize_pops, POP2 = symp_mexicana_pops), # mexicana ancestry within mexicana (outliers defined by local maize)
        expand("diversity/results/pi/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/mexicana/{POP2}.4pop.outliers{POP1}.pi.allChr.pestPG", zip, POP1 = symp_maize_pops, POP2 = symp_mexicana_pops),
        # fst between mexicana ancestry within sympatric mexicana and within local sympatric maize (at the introgression outliers for the local maize)
        expand("diversity/results/fst/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/mexicana/{POP1}.{POP2}.1pop.outliers{POP1}.fst.allChr.txt", zip, POP1 = symp_maize_pops, POP2 = symp_mexicana_pops),
        expand("diversity/results/fst/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/mexicana/{POP1}.{POP2}.4pop.outliers{POP1}.fst.allChr.txt", zip, POP1 = symp_maize_pops, POP2 = symp_mexicana_pops),
        expand("diversity/results/fst/" + prefix_all + "/Ne{Ne}_{YESNO}Boot/HOMOZYG/summary_pop_pairs_fst.outliers.allChr.txt", Ne = 10000, YESNO = "yes")
    params:
        p = "med2"
    resources:
        time_min = 30,
        mem = 2

rule fst_maize_anc_outliers:
    input:
        expand("diversity/results/pi/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/maize/{POP}.{n}pop.outliers{POP}.pi.allChr.pestPG", POP = symp_mexicana_pops, n = [1, 4]), # introgressed maize tracts within sympatric mexicana
        expand("diversity/results/pi/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/maize/{POP2}.1pop.outliers{POP1}.pi.allChr.pestPG", zip, POP1 = symp_mexicana_pops, POP2 = symp_maize_pops), # maize ancestry within maize (outliers defined by local maize)
        expand("diversity/results/pi/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/maize/{POP2}.4pop.outliers{POP1}.pi.allChr.pestPG", zip, POP1 = symp_mexicana_pops, POP2 = symp_maize_pops),
        # fst between maize ancestry within sympatric mexicana and within local sympatric maize (at the introgression outliers for the local mexicana)
        expand("diversity/results/fst/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/maize/{POP1}.{POP2}.1pop.outliers{POP1}.fst.allChr.txt", zip, POP1 = symp_mexicana_pops, POP2 = symp_maize_pops),
        expand("diversity/results/fst/" + prefix_all + "/Ne10000_yesBoot/HOMOZYG/maize/{POP1}.{POP2}.4pop.outliers{POP1}.fst.allChr.txt", zip, POP1 = symp_mexicana_pops, POP2 = symp_maize_pops)
    params:
        p = "med2"
    resources:
        time_min = 30,
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
