import os
import numpy as np
from csv import DictReader
import itertools

#configfile: "config.yaml"
path_hilo = os.getcwd() + "/" # get absolute path to this hilo git directory on the local machine

# wildcards
wildcard_constraints:
    ID = "[A-Za-z0-9]+",
    POP = "pop[0-9]+|allopatric_maize|allopatric_maize_subsample[0-9]+|sympatric_maize|sympatric_mexicana|allopatric_mexicana|parv", # for f4 stats allopatric_maize needs to be treated like a pop
    NOTPOP = "not[0-9]+", # eg. not360 is all maize EXCEPT population 360
    POP1 = "pop[0-9]+|not[0-9]+", # pop 1 and pop 2 for pairwise Fst calculations and saf/sfs. Can be 1 population (e.g. pop360) or all individuals in that subspecies except 1 population (e.g. not360)
    POP2 = "pop[0-9]+|not[0-9]+",
    BOOT = "[0-9]+", # bootstrap
    YESNO = "yes|no",
    Ne = "[0-9]+",
    POSNEG = "pos|neg",
    SIG = "fdr05|perc05|p05", # outlier significance cutoffs: 5% FDR, 5% empirical cutoff, p = 5% (5% percentile from null simulations)
    STAT = "meanAnc|lmElev|maize_anc|mexicana_anc|parv_anc", # statistics defining outliers
    GROUP = "sympatric_maize|sympatric_mexicana|allopatric_maize|allopatric_mexicana|parv",
    ZEA = "maize|mexicana|parv",
    ANCESTRY = "maize|mexicana|parv",
    ALLO_MEX = "allopatric_maize|pop22", # used in f4s. pop22 is Amecameca (at 2467m, the highest elevation of 3 allopatric pops)
    SUBSAMPLE = "[0-9]+",
    REGION = "region_[0-9]+",
    r = "[0-9]+", # numeric only
    FEATURE = "r|cd|frac", # recombination rate cM/Mb (r), coding bp/cM (cd), or frac coding bp (frac)
    COVERAGE = "ALL|Over0.5x", # we only estimate local ancestry for individuals with >0.5x mean coverage. whereas allopatric pops and global ancestry estimates we used a less stringent cutoff (individuals with >0.1x coverage)
    WIN = "[0-9]+", # window size
    STEP = "[0-9]+", # window step size (non-overlapping windows if STEP = WIN)
    n = "[0-9]+",
    PREFIX = "HILO_MAIZE55|HILO_MAIZE55_PARV50",
    K = "[0-9]+"

# reference genome and associated files
#ref = "/home/ecalfee/hilo/data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
ref = path_hilo + "data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
fai = ref + ".fai"
ref_chr = path_hilo + "data/refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
ref_anno = path_hilo + "data/refMaize/geneAnnotations/Zea_mays.B73_RefGen_v4.41.chr.gff3"
trip_anc = "filtered_bams/results/SRR7758238/TRIP.fa.gz"

# recombination map
rmap_ext = path_hilo + "linkage_map/results/ogut_2015_rmap_v2_to_v4_EXTENDED.txt"

# 1cM genomic windows for whole genome (1477 windows)
windows_1cM = []
#with open("ancestry_by_r/results/map_pos_1cM_windows.txt", "r") as read_obj:
#    csv_dict_reader = DictReader(read_obj, delimiter = "\t")
#    for row in csv_dict_reader:
#        windows_1cM.append(row["window"])
for i in range(1, 1478):
    windows_1cM.append("W" + str(i))


# all bams included in combined sample (est_coverage > 0.05x)
prefix_all = "HILO_MAIZE55"
#prefix_all = "TEST"

# list of all included bams and ids (over minimum 0.05x coverage)
with open("samples/HILO_MAIZE55_bams.list") as f:
    mex_maize_bams = f.read().splitlines()
with open("samples/HILO_MAIZE55_PARV50_bams.list") as f:
    mex_maize_parv_bams = f.read().splitlines()
with open("samples/HILO_MAIZE55_ids.list") as f:
    mex_maize_ids = f.read().splitlines()
with open("samples/HILO_MAIZE55_PARV50_ids.list") as f:
    mex_maize_parv_ids = f.read().splitlines()
# add bams for outgroup tripsacum
with open("samples/ALL_byPop/trip_bams.list") as f:
    trip_bams = f.read().splitlines()
# add bams for parviglumis
with open("samples/ALL_byPop/parv_bams.list") as f:
    parv_bams = f.read().splitlines()
all_bams_mex_maize_parv_trip = mex_maize_parv_bams + trip_bams

# define functions to get list of input bams and bais
def get_all_bams(prefix):
    with open(path_hilo + "samples/" + prefix + "_bams.list") as f:
        bams = f.read().splitlines()
    return bams

def get_all_bais(prefix):
    with open(path_hilo + "samples/" + prefix + "_bams.list") as f:
        bams = f.read().splitlines()
        bais = [bam + ".bai" for bam in bams]
    return bais

def get_all_ids(prefix):
    with open(path_hilo + "samples/" + prefix + "_ids.list") as f:
        ids = f.read().splitlines()
    return ids

# function to look up the path to a single bam using the prefix for the set of bams and the sample id
def lookup_one_bam(prefix, id):
    all_ids = get_all_ids(prefix)
    all_bams = get_all_bams(prefix)
    my_bam = all_bams[all_ids.index(id)]
    return my_bam

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
#include: "map/Snakefile"
#include: "filtered_bams/Snakefile"
#include: "variant_sites/Snakefile"
#include: "global_ancestry/Snakefile"
#include: "linkage_map/Snakefile"
#include: "local_ancestry/Snakefile"
include: "ancestry_by_r/Snakefile"
include: "ZAnc/Snakefile"
include: "diversity/Snakefile"
include: "mhl1_inv/Snakefile"
#include: "domestication_scan/Snakefile"
#include: "wavelets/Snakefile"

## all:  main rule to run all workflows: main figures, tables supplement, figures supplement, other output
rule all:
    input:
        # main figures:
        # Fig 1
        "map/plots/mexico_lines_elev_teo_color.png",
        "../hilo_manuscript/figures_main/mexico_lines_elev_teo_color.tif",

        # Fig 2
        #"global_ancestry/plots/global_anc_multi.png",
        #"../hilo_manuscript/figures_main/global_anc_multi.tif",
        "../hilo/global_ancestry/plots/HILO_MAIZE55_PARV50_global_anc_multi_K3.png",
        "../hilo_manuscript/figures_main/HILO_MAIZE55_PARV50_global_anc_multi_K3.tif",

        # Fig 3
        #"diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/fst_within_maize_or_mexicana_ancestry_genomewide_heatmap_both.png",
        #"../hilo_manuscript/figures_main/Ne10000_yesBoot_fst_within_maize_or_mexicana_ancestry_genomewide_heatmap_both.tif",
        "../hilo_manuscript/figures_main/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_fst_within_maize_or_mexicana_ancestry_genomewide_heatmap_both.tif",
        "diversity/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_fst_within_maize_or_mexicana_ancestry_genomewide_heatmap_both.png",

        # Fig 4
        #"ancestry_by_r/plots/K2_by_r_multi_panel.png",
        #"../hilo_manuscript/figures_main/K2_by_r_multi_panel.tif",
        "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_by_r_multi_panel.png",
        "../hilo_manuscript/figures_main/HILO_MAIZE55_PARV50_K3_by_r_multi_panel.tif",

        # Fig 5
        #"ZAnc/plots/Ne10000_yesBoot/maize_shared_outliers_chr_4.png",
        #"../hilo_manuscript/figures_main/Ne10000_yesBoot_maize_shared_outliers_chr_4.tif",
        "ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_maize_shared_outliers_chr_4.png",
        "../hilo_manuscript/figures_main/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_maize_shared_outliers_chr_4.tif",

        # Fig 6
        #"ZAnc/plots/Ne10000_yesBoot/network_peak_sharing_data_only.png",
        #"../hilo_manuscript/figures_main/Ne10000_yesBoot_network_peak_sharing_data_only.tif",
        "ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_network_peak_sharing_data_only.png",
        "../hilo_manuscript/figures_main/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_network_peak_sharing_data_only.tif",

        # Fig 7
        #"ZAnc/plots/Ne10000_yesBoot/multi_maize_mexicana_genome_scan.png",
        #"../hilo_manuscript/figures_main/Ne10000_yesBoot_multi_maize_mexicana_genome_scan.tif",
        "ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_multi_maize_mexicana_genome_scan.png",
        "../hilo_manuscript/figures_main/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_multi_maize_mexicana_genome_scan.tif",


        # tables supplement
        #samples/population_metadata.csv, # made outside snakemake pipeline (also saved to github)
        #samples/parviglumis_50_SRA_IDs.csv, # made outside snakemake pipeline (also saved to github)

        #"ancestry_by_r/tables/spearmans_rho_ngsadmix.tex",
        #"../hilo_manuscript/tables/spearmans_rho_ngsadmix.tex",
        "ancestry_by_r/tables/HILO_MAIZE55_PARV50_K3_spearmans_rho_ngsadmix.tex",
        "../hilo_manuscript/tables/HILO_MAIZE55_PARV50_K3_spearmans_rho_ngsadmix.tex",

        "ancestry_by_r/tables/spearmans_rho_f4_sympatric_maize_pop22.tex",
        "../hilo_manuscript/tables/spearmans_rho_f4_sympatric_maize_pop22.tex",

        #"ancestry_by_r/tables/elev_r_interaction_5.tex",
        #"../hilo_manuscript/tables/elev_r_interaction_5.tex",
        "ancestry_by_r/tables/HILO_MAIZE55_PARV50_K3_elev_r_interaction_5.tex",
        "../hilo_manuscript/tables/HILO_MAIZE55_PARV50_K3_elev_r_interaction_5.tex",

        "ancestry_by_r/tables/spearmans_rho_f4_sympatric_mexicana_pop22.tex",
        "../hilo_manuscript/tables/spearmans_rho_f4_sympatric_mexicana_pop22.tex",

        #"ancestry_by_r/tables/spearmans_rho_local_ancestry.tex",
        #"../hilo_manuscript/tables/spearmans_rho_local_ancestry.tex",
        "ancestry_by_r/tables/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_spearmans_rho_local_ancestry.tex",
        "../hilo_manuscript/tables/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_spearmans_rho_local_ancestry.tex",

        #"domestication_scan/tables/HILO_MAIZE55/Ne10000_yesBoot/domestication_genes.tex",
        #"../hilo_manuscript/tables/Ne10000_yesBoot_domestication_genes.tex",
        "ZAnc/tables/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_genes_mapped_to_outliers.tex",
        "../hilo_manuscript/tables/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_genes_mapped_to_outliers.tex",


        # figures supplement
        #"global_ancestry/plots/HILO_MAIZE55_pca.png",
        #"../hilo_manuscript/figures_supp/HILO_MAIZE55_pca.tif",
        "global_ancestry/plots/HILO_MAIZE55_PARV50_pca.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_pca.tif",

        # K=3 structure results as a barplot for each individual
        "global_ancestry/plots/HILO_MAIZE55_PARV50_structure_K3.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_structure_K3.tif",

        # supplementary plot of within-parviglumis-ancestry diversity (fst) as a heat map
        "diversity/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_fst_within_parviglumis_ancestry_genomewide_heatmap_both.png",
        "../hilo_manuscript/figures_main/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_fst_within_parviglumis_ancestry_genomewide_heatmap_both.tif",

        #"local_ancestry/plots/admix_times_Ne10000_yesBoot.png",
        #"../hilo_manuscript/figures_supp/admix_times_Ne10000_yesBoot.tif",
        "local_ancestry/plots/HILO_MAIZE55_PARV50_K3_admix_times_Ne10000_yesBoot.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_admix_times_Ne10000_yesBoot.tif",

        #"diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/pi_within_mexicana_ancestry_peaks.png",
        #"../hilo_manuscript/figures_supp/Ne10000_yesBoot_pi_within_mexicana_ancestry_peaks.tif",
        "diversity/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_pi_within_mexicana_ancestry_peaks.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_pi_within_mexicana_ancestry_peaks.tif",

        #"diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/pi_within_maize_ancestry.png",
        #"../hilo_manuscript/figures_supp/Ne10000_yesBoot_pi_within_maize_ancestry.tif",
        "diversity/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_pi_within_maize_ancestry.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_pi_within_maize_ancestry.tif",

        # tree_f4_stats.png, # made outside of snakemake pipeline

        "ancestry_by_r/plots/f4_sympatric_maize_pop22_byr5.png",
        "../hilo_manuscript/figures_supp/f4_sympatric_maize_pop22_byr5.tif",

        "ancestry_by_r/plots/f4_sympatric_maize_pop22_bycd5.png",
        "../hilo_manuscript/figures_supp/f4_sympatric_maize_pop22_bycd5.tif",

        #"ancestry_by_r/plots/K2_by_r_bootstrap_lm_elevation_color_elev.png",
        #"../hilo_manuscript/figures_supp/K2_by_r_bootstrap_lm_elevation_color_elev.tif",
        "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_by_r_bootstrap_lm_elevation_color_elev.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_by_r_bootstrap_lm_elevation_color_elev.tif",

        # supplementary plot of all NGSAdmix-estimated ancestries, including parviglumis ancestry, by recombination rate
        "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_by_r_bootstrap_sympatric_and_allopatric.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_by_r_bootstrap_sympatric_and_allopatric.tif",

        #"ancestry_by_r/plots/K2_by_cd_bootstrap_sympatric_and_allopatric.png",
        #"../hilo_manuscript/figures_supp/K2_by_cd_bootstrap_sympatric_and_allopatric.tif",
        "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_by_cd_bootstrap_sympatric_and_allopatric.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_by_cd_bootstrap_sympatric_and_allopatric.tif",

        "ancestry_by_r/plots/f4_sympatric_mexicana_pop22_byr5.png",
        "../hilo_manuscript/figures_supp/f4_sympatric_mexicana_pop22_byr5.tif",

        "ancestry_by_r/plots/f4_sympatric_mexicana_pop22_bycd5.png",
        "../hilo_manuscript/figures_supp/f4_sympatric_mexicana_pop22_bycd5.tif",

        #"ancestry_by_r/plots/local_anc_by_r_continuous.png",
        #"../hilo_manuscript/figures_supp/local_anc_by_r_continuous.tif",
        "ancestry_by_r/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_local_anc_by_r_continuous.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_local_anc_by_r_continuous.tif",

        #"ZAnc/plots/Ne10000_yesBoot/combmatrix_peak_sharing_maize.png",
        #"../hilo_manuscript/figures_supp/Ne10000_yesBoot_combmatrix_peak_sharing_maize.tif",
        "ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_combmatrix_peak_sharing_maize.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_combmatrix_peak_sharing_maize.tif",

        #"ZAnc/plots/Ne10000_yesBoot/combmatrix_peak_sharing_mexicana.png",
        #"../hilo_manuscript/figures_supp/Ne10000_yesBoot_combmatrix_peak_sharing_mexicana.tif",
        "ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_combmatrix_peak_sharing_mexicana.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_combmatrix_peak_sharing_mexicana.tif",

        # local ancestry across all the individual chromosomes for maize and mexicana
        #expand("ZAnc/plots/Ne10000_yesBoot/{ZEA}_shared_outliers_chr_{i}.png", i = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], ZEA = zea),
        #expand("../hilo_manuscript/figures_supp/Ne10000_yesBoot_{ZEA}_shared_outliers_chr_{i}.tif", i = [1, 2, 3, 5, 6, 7, 8, 9, 10], ZEA = "maize"), # skips chr4 (that's a main figure)
        #expand("../hilo_manuscript/figures_supp/Ne10000_yesBoot_{ZEA}_shared_outliers_chr_{i}.tif", i = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], ZEA = "mexicana"),
        expand("ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_{ZEA}_shared_outliers_chr_{i}.png", i = [1, 2, 3, 5, 6, 7, 8, 9, 10], ZEA = "maize"),
        expand("ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_{ZEA}_shared_outliers_chr_{i}.png", i = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], ZEA = "mexicana"),
        expand("../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_{ZEA}_shared_outliers_chr_{i}.tif", i = [1, 2, 3, 5, 6, 7, 8, 9, 10], ZEA = "maize"),
        expand("../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_{ZEA}_shared_outliers_chr_{i}.tif", i = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], ZEA = "mexicana"),

        #"diversity/plots/HILO_MAIZE55/Ne10000_yesBoot/local_fst_within_mexicana_ancestry_peaks.png",
        #"../hilo_manuscript/figures_supp/Ne10000_yesBoot_local_fst_within_mexicana_ancestry_peaks.tif",
        "diversity/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_local_fst_within_mexicana_ancestry_peaks.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_local_fst_within_mexicana_ancestry_peaks.tif",

        #"ZAnc/plots/Ne10000_yesBoot/QQ.png",
        #"../hilo_manuscript/figures_supp/Ne10000_yesBoot_QQ.tif",
        "ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_QQ.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_QQ.tif",

        # supplementary plots of parviglumis ancestry introgression peaks
        "ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_supp_parv_genome_scan.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_supp_parv_genome_scan.tif",

        #"mhl1_inv/plots/HILO_MAIZE55/K2/Ne10000_yesBoot/mhl1_inv_ancestry.png",
        #"../hilo_manuscript/figures_supp/HILO_MAIZE55_K2_Ne10000_yesBoot_mhl1_inv_ancestry.tif",
        "mhl1_inv/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_mhl1_inv_ancestry.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_mhl1_inv_ancestry.tif",

        #"mhl1_inv/plots/HILO_MAIZE55/K2/Ne10000_yesBoot/mhl1_inv_pca.png",
        #"../hilo_manuscript/figures_supp/HILO_MAIZE55_K2_Ne10000_yesBoot_mhl1_inv_pca.tif",
        "mhl1_inv/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_mhl1_inv_pca.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_mhl1_inv_pca.tif",

        # supplementary plot of hpc1
        "ZAnc/plots/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_hpc1_lm_ancestry.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_Ne10000_yesBoot_hpc1_lm_ancestry.tif",

        "filtered_bams/plots/p_seq_counts.png",
        "../hilo_manuscript/figures_supp/p_seq_counts.tif",

        "linkage_map/plots/ogut_2015_v2_to_v4_rmap.png",
        "../hilo_manuscript/figures_supp/ogut_2015_v2_to_v4_rmap.tif",

        # Sensitivity to Ne
        #expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55_PARV50/K3/Ne{Ne}_noBoot/anc/{ZEA}.{ANCESTRY}_anc.bed", Ne = Nes, ZEA = zea, ANCESTRY = ["mexicana", "maize", "parv"]),
        "local_ancestry/plots/HILO_MAIZE55_PARV50_K3_sensitivity_to_Ne_admix_times.png",
        "../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_sensitivity_to_Ne_admix_times.tif",
        "local_ancestry/results/HILO_MAIZE55_PARV50_K3_sensitivity_to_Ne_local_ancestry.txt",
        "../hilo_manuscript/tables/HILO_MAIZE55_PARV50_K3_sensitivity_to_Ne_local_ancestry.tex",

        # intermediate quality control reports:
        "filtered_bams/metrics/fastQC/multiqc/multiqc_report.html", # bam metrics files
        "filtered_bams/metrics/fastQC_trimmed/multiqc/multiqc_report.html",
        "filtered_bams/metrics/picard/multiqc/multiqc_report.html",
        "filtered_bams/metrics/flagstat/multiqc/multiqc_report.html",

        # other output files/results that go into the manuscript text:
        # lm results for mexicana ancestry globally ~ elevation:
        "global_ancestry/tables/HILO_MAIZE55_PARV50_lm_elevation_K3.tex",

        # summary of admixture times (median, min, max) estimated by ancestry_hmm:
        "local_ancestry/results/HILO_MAIZE55_PARV50_K3_admix_times_Ne10000_yesBoot.txt",

        # stats for percent of MVN simulations truncated to be within [0,1] range
        expand("ZAnc/results/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/{ZEA}.MVN.truncated.stats.txt", ZEA = zea),

        # stats for number of bootstraps with ambiguous K ancestry assignments (ancestry ~ r/cd analysis)
        expand("ancestry_by_r/results/bootstrap_1cM/HILO_MAIZE55_PARV50/{FEATURE}5_K3_boot_drop_stats.txt", FEATURE = ["r", "cd"]),

        # outliers for maize, mexicana or parviglumis outliers as bed file
        expand("ZAnc/results/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/{ZEA}_pos_{ANCESTRY}_anc_outliers.fdr05.bed", ZEA = zea, ANCESTRY = ["mexicana", "maize", "parv"]),

        # outliers for lm mexicana ancestry ~ elevation as bed file
        expand("ZAnc/results/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/{ZEA}_pos_lmElev_outliers.fdr05.bed", ZEA = zea),

        # count of flowering time related genes overlapping outliers from genome scan
        #"ZAnc/results/HILO_MAIZE55/Ne10000_yesBoot/flowering_time_genes_v4.plus20kb.overlap.summary_overlap_outliers.txt"
        "ZAnc/results/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/flowering_time_genes_v4.plus20kb.overlap.summary_overlap_outliers.txt"

    params:
        p = "med2"
    resources:
        time_min = 60,
        mem = 8

## some: for running a subset of analyses
rule some:
    input:
        expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55_PARV50/K3/Ne{Ne}_noBoot/anc/{ZEA}/{POP}.anc.freq", Ne = Nes, POP = symp_pops, ZEA = ["maize", "mexicana", "parv"]),
        expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/K2/Ne{Ne}_noBoot/anc/{ZEA}/{POP}.anc.freq", Ne = Nes, POP = symp_pops, ZEA = zea),

        expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55_PARV50/K3/Ne{Ne}_noBoot/anc/{ZEA}/{POP}.alpha.ind", Ne = Nes, POP = symp_pops, ZEA = ["maize", "mexicana", "parv"]),
        expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55_PARV50/K3/Ne{Ne}_noBoot/anc/{ZEA}/{POP}.alpha.ind", Ne = Nes, POP = symp_pops, ZEA = ["maize", "mexicana"]),

        expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55_PARV50/K3/Ne{Ne}_yesBoot/anc/{ZEA}/{POP}.alpha.ind", Ne = 10000, POP = symp_pops, ZEA = ["maize", "mexicana", "parv"]),
        expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/K2/Ne{Ne}_yesBoot/anc/{ZEA}/{POP}.alpha.ind", Ne = 10000, POP = symp_pops, ZEA = ["maize", "mexicana"]),

        #expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55_PARV50/K3/Ne{Ne}_yesBoot/anc/{ZEA}/{POP}.anc.freq", Ne = 10000, POP = symp_pops, ZEA = ["maize", "mexicana", "parv"]), # only bootstrap t for Ne=10000
        #expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/K2/Ne{Ne}_yesBoot/anc/{ZEA}/{POP}.anc.freq", Ne = Nes, POP = symp_pops, ZEA = zea),

        #expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/anc/{ZEA}.pops.anc.RData", ZEA = zea),
        #expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/K2/Ne10000_yesBoot/anc/{ZEA}.pops.anc.RData", ZEA = zea),

        expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55_PARV50/K3/Ne{Ne}_noBoot/{POP}.times", Ne = Nes, POP = symp_pops),
        expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/K2/Ne{Ne}_noBoot/{POP}.times", Ne = Nes, POP = symp_pops),
        expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55_PARV50/K3/Ne{Ne}_yesBoot/{POP}.times", Ne = 10000, POP = symp_pops),
        expand("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/K2/Ne{Ne}_yesBoot/{POP}.times", Ne = 10000, POP = symp_pops),
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_admix_times_Ne10000_yesBoot.tif", zip, PREFIX = ["HILO_MAIZE55", "HILO_MAIZE55_PARV50"], K = [2, 3]),

        #expand("ZAnc/results/{PREFIX}/K{K}/Ne10000_yesBoot/{ZEA}.MVN.RData", PREFIX = "HILO_MAIZE55", K = 2, ZEA = zea),
        #expand("ZAnc/results/{PREFIX}/K{K}/Ne10000_yesBoot/{ZEA}.MVN.RData", PREFIX = "HILO_MAIZE55_PARV50", K = 3, ZEA = zea),

        expand("ZAnc/results/{PREFIX}/K{K}/Ne10000_yesBoot/{ZEA}.MVN.truncated.stats.txt", PREFIX = "HILO_MAIZE55", K = 2, ZEA = zea),
        expand("ZAnc/results/{PREFIX}/K{K}/Ne10000_yesBoot/{ZEA}.MVN.truncated.stats.txt", PREFIX = "HILO_MAIZE55_PARV50", K = 3, ZEA = zea),

        #expand("ZAnc/results/{PREFIX}/K{K}/Ne10000_yesBoot/{ZEA}.lmElev.fit.RData", PREFIX = "HILO_MAIZE55", K = 2, ZEA = zea),
        #expand("ZAnc/results/{PREFIX}/K{K}/Ne10000_yesBoot/{ZEA}.lmElev.fit.RData", PREFIX = "HILO_MAIZE55_PARV50", K = 3, ZEA = zea),
        #expand("ZAnc/results/{PREFIX}/K{K}/Ne10000_yesBoot/{ZEA}.meanAnc.fdr.RData", PREFIX = "HILO_MAIZE55", K = 2, ZEA = zea),
        #expand("ZAnc/results/{PREFIX}/K{K}/Ne10000_yesBoot/{ZEA}.meanAnc.fdr.RData", PREFIX = "HILO_MAIZE55_PARV50", K = 3, ZEA = zea),

        expand("ZAnc/results/{PREFIX}/K{K}/Ne10000_yesBoot/{ZEA}_pos_{ANCESTRY}_anc_outliers.fdr05.bed", PREFIX = "HILO_MAIZE55", K = 2, ZEA = zea, ANCESTRY = zea),
        expand("ZAnc/results/{PREFIX}/K{K}/Ne10000_yesBoot/{ZEA}_pos_{ANCESTRY}_anc_outliers.fdr05.bed", PREFIX = "HILO_MAIZE55_PARV50", K = 3, ZEA = zea, ANCESTRY = ["mexicana", "maize", "parv"]),
        expand("ZAnc/results/{PREFIX}/K{K}/Ne10000_yesBoot/{ZEA}_pos_lmElev_outliers.fdr05.bed", PREFIX = "HILO_MAIZE55", K = 2, ZEA = zea),
        expand("ZAnc/results/{PREFIX}/K{K}/Ne10000_yesBoot/{ZEA}_pos_lmElev_outliers.fdr05.bed", PREFIX = "HILO_MAIZE55", K = 3, ZEA = "mexicana"),

        expand("ZAnc/results/{PREFIX}/K{K}/Ne10000_yesBoot/flowering_time_genes_v4.plus20kb.overlap.summary_overlap_outliers.txt", zip, PREFIX = ["HILO_MAIZE55", "HILO_MAIZE55_PARV50"], K = [2, 3]),

        expand("ZAnc/plots/{PREFIX}_K{K}_Ne10000_yesBoot_{ZEA}_mean_{ANCESTRY}_anc.png", PREFIX = "HILO_MAIZE55_PARV50", K = 3, ZEA = zea, ANCESTRY = ["mexicana", "maize", "parv"]),
        expand("ZAnc/plots/{PREFIX}_K{K}_Ne10000_yesBoot_{ZEA}_mean_{ANCESTRY}_anc.png", PREFIX = "HILO_MAIZE55", K = 2, ZEA = zea, ANCESTRY = ["mexicana", "maize"]),
        expand("ZAnc/plots/{PREFIX}_K{K}_Ne10000_yesBoot_{ZEA}_slope_elev.png", PREFIX = "HILO_MAIZE55_PARV50", K = 3, ZEA = zea),
        expand("ZAnc/plots/{PREFIX}_K{K}_Ne10000_yesBoot_{ZEA}_slope_elev.png", PREFIX = "HILO_MAIZE55", K = 2, ZEA = zea),
        expand("ZAnc/plots/{PREFIX}_K{K}_Ne10000_yesBoot_multi_maize_mexicana_genome_scan.png", zip, PREFIX = ["HILO_MAIZE55", "HILO_MAIZE55_PARV50"], K = [2, 3]),

        expand("../hilo_manuscript/figures_supp/{PREFIX}_lm_{ANCESTRY}_by_pop_elevation_K3.tif", PREFIX = "HILO_MAIZE55_PARV50", ANCESTRY = ["maize", "parviglumis"]),
        "global_ancestry/tables/HILO_MAIZE55_PARV50_lm_elevation_K3.tex",


        expand("ancestry_by_r/results/bootstrap_1cM/{PREFIX}/{FEATURE}5_K{K}.Rdata", PREFIX = "HILO_MAIZE55", K = 2, FEATURE = ["cd", "r"]),
        expand("ancestry_by_r/results/bootstrap_1cM/{PREFIX}/{FEATURE}5_K{K}.Rdata", PREFIX = "HILO_MAIZE55_PARV50", K = 3, FEATURE = ["cd", "r"]),
        expand("ancestry_by_r/results/local_anc_1cM/{PREFIX}/K{K}/Ne10000_yesBoot/{ANCESTRY}_anc/{POP}.anc.wind", PREFIX = "HILO_MAIZE55_PARV50", K = 3, ANCESTRY = ["mexicana", "maize", "parv"], POP = symp_pops),
        expand("ancestry_by_r/results/local_anc_1cM/{PREFIX}/K{K}/Ne10000_yesBoot/{ANCESTRY}_anc/{POP}.anc.wind", PREFIX = "HILO_MAIZE55", K = 2, ANCESTRY = ["mexicana", "maize"], POP = symp_pops),

        # ancestry ~ r (or coding density) plots and tables: NGSAdmix results
        expand("../hilo_manuscript/figures_main/{PREFIX}_K{K}_by_r_multi_panel.tif", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_by_r_bootstrap_sympatric_and_allopatric.tif", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_by_cd_bootstrap_sympatric_and_allopatric.tif", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_by_r_bootstrap_lm_elevation_color_elev.tif", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        expand("../hilo_manuscript/tables/{PREFIX}_K{K}_elev_r_interaction_5.tex", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        expand("../hilo_manuscript/tables/{PREFIX}_K{K}_spearmans_rho_ngsadmix.tex", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        # ancestry ~ r (or coding density) plots and tables: local ancestry results
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_Ne10000_yesBoot_local_anc_by_r_continuous.tif", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        expand("ancestry_by_r/tables/{PREFIX}_K{K}_Ne10000_yesBoot_spearmans_rho_local_ancestry.tex", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        expand("../hilo_manuscript/tables/{PREFIX}_K{K}_Ne10000_yesBoot_spearmans_rho_local_ancestry.tex", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        # ancestry ~ r (or coding density) plots and tables: f4 results (not changed by K=3 choice)

        expand("mhl1_inv/plots/{PREFIX}_K{K}_Ne10000_yesBoot_mhl1_inv_pca.png", zip, PREFIX = ["HILO_MAIZE55_PARV50", "HILO_MAIZE55"], K = [3,2]),
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_Ne10000_yesBoot_mhl1_inv_pca.tif", zip, PREFIX = ["HILO_MAIZE55_PARV50", "HILO_MAIZE55"], K = [3,2]),
        expand("mhl1_inv/plots/{PREFIX}_K{K}_Ne10000_yesBoot_mhl1_inv_ancestry.png", zip, PREFIX = ["HILO_MAIZE55_PARV50", "HILO_MAIZE55"], K = [3,2]),
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_Ne10000_yesBoot_mhl1_inv_ancestry.tif", zip, PREFIX = ["HILO_MAIZE55_PARV50", "HILO_MAIZE55"], K = [3,2]),

        #"../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_sensitivity_to_Ne_admix_times.tif",
        #"../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_sensitivity_to_Ne_admix_times.tif",
        #"local_ancestry/plots/HILO_MAIZE55_PARV50_K3_sensitivity_to_Ne_local_ancestry.png",
        #"../hilo_manuscript/figures_supp/HILO_MAIZE55_PARV50_K3_sensitivity_to_Ne_local_ancestry.tif",

        # introgression peaks
        expand("../hilo_manuscript/figures_main/{PREFIX}_K{K}_Ne10000_yesBoot_network_peak_sharing_data_only.tif", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_Ne10000_yesBoot_combmatrix_peak_sharing_maize.tif", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_Ne10000_yesBoot_combmatrix_peak_sharing_mexicana.tif", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        expand("ZAnc/plots/{PREFIX}_K{K}_Ne10000_yesBoot/{ZEA}_shared_outliers_chr_{i}.png", i = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], ZEA = zea, PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        expand("../hilo_manuscript/figures_main/{PREFIX}_K{K}_Ne10000_yesBoot_maize_shared_outliers_chr_4.tif", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_Ne10000_yesBoot_{ZEA}_shared_outliers_chr_{i}.tif", i = [1, 2, 3, 5, 6, 7, 8, 9, 10], ZEA = "maize", PREFIX = "HILO_MAIZE55_PARV50", K = 3),
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_Ne10000_yesBoot_{ZEA}_shared_outliers_chr_{i}.tif", i = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], ZEA = "mexicana", PREFIX = "HILO_MAIZE55_PARV50", K = 3)

    params:
        p = "med2"
    resources:
        time_min = 60,
        mem = 8

rule files_for_jeff:
    input:
        expand("wavelets/results/alleleFreqs/HILO_MAIZE55_PARV50/K3/{POP}.mafs.gz", POP = symp_pops + allo_mex_pops + ["allopatric_maize", "allopatric_mexicana", "parv"]),
        expand("wavelets/results/alleleFreqs/HILO_MAIZE55/K2/{POP}.mafs.gz", POP = symp_pops + allo_mex_pops + ["allopatric_maize", "allopatric_mexicana"])
    params:
        p = "med2"
    resources:
        time_min = 60,
        mem = 2

rule diversity:
    input:
        "diversity/results/fst/HILO_MAIZE55/K2/Ne10000_yesBoot/HOMOZYG/summary_pop_pairs_fst.allChr.txt",
        "diversity/results/fst/HILO_MAIZE55_PARV50/K3/Ne10000_yesBoot/HOMOZYG/summary_pop_pairs_fst.allChr.txt",
        expand("../hilo_manuscript/figures_main/{PREFIX}_K{K}_Ne10000_yesBoot_fst_within_maize_or_mexicana_ancestry_genomewide_heatmap_both.tif", zip, PREFIX = ["HILO_MAIZE55", "HILO_MAIZE55_PARV50"], K = [2, 3]),
        expand("diversity/plots/{PREFIX}_K{K}_Ne10000_yesBoot_fst_within_maize_or_mexicana_ancestry_genomewide_heatmap_both.png", zip, PREFIX = ["HILO_MAIZE55", "HILO_MAIZE55_PARV50"], K = [2, 3]),
        expand("diversity/results/pi/{PREFIX}/K{K}/Ne{Ne}_{YESNO}Boot/HOMOZYG/summary_pop_pi.allChr.txt", PREFIX = "HILO_MAIZE55_PARV50", K = 3, Ne = 10000, YESNO = "yes"),
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_Ne{Ne}_{YESNO}Boot_pi_within_mexicana_ancestry_peaks.tif", PREFIX = "HILO_MAIZE55_PARV50", K = 3, Ne = 10000, YESNO = "yes"),
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_Ne{Ne}_{YESNO}Boot_pi_within_maize_ancestry.tif", PREFIX = "HILO_MAIZE55_PARV50", K = 3, Ne = 10000, YESNO = "yes"),
        expand("../hilo_manuscript/figures_supp/{PREFIX}_K{K}_Ne{Ne}_{YESNO}Boot_local_fst_within_mexicana_ancestry_peaks.tif", PREFIX = "HILO_MAIZE55_PARV50", K = 3, Ne = 10000, YESNO = "yes")
    params:
        p = "med2"
    resources:
        time_min = 60,
        mem = 2

rule parv:
    input:
        expand("local_ancestry/results/alloFreqs/{PREFIX}/{GROUP}/{REGION}.mafs.gz", PREFIX = "HILO_MAIZE55_PARV50", GROUP = "parv", REGION = list(regions_dict.keys()))
    params:
        p = "med2"
    resources:
        time_min = 60,
        mem = 2
