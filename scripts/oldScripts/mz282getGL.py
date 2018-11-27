#!/usr/bin/python
# script to parse maize 282 panel and extract genotype likelihoods
# this script should be used in the python virtual environment jPythonEnv
# which can be activated by jPythonEnv$ source bin/activate
# and deactivated using anyDirectory$ deactivate
# RUN on Python 3
import sys
import vcf
import csv
import math


# helper function to take genotype likelihood 3 fields
# and convert from normalized & thread scaled to 
# on in the case that GL is just None (no genotype) .. outputs (.333333, .333333, .333333) for new beagle-format GL
def convertGL(phredNLL): # takes in a phred normalized genotype likelihood = -10*log(lik(data|genotype)/lik(data|most likely genotype))
    if phredNLL is None: # not sequenced at this base -- every genotype equally likely to produce null data
        GL = [1.0/3.0, 1.0/3.0, 1.0/3.0]
    else:
        A = math.exp(-0.1*phredNLL[0]) # GL for genotype 0
        B = math.exp(-0.1*phredNLL[1]) # GL for genotype 1
        C = math.exp(-0.1*phredNLL[2]) # GL for genotype 2
        # normalize so they sum to one
        GL = [A/(A+B+C), B/(A+B+C), C/(A+B+C)]
    GL = ["%.6f"%f for f in GL] # round to 6 decimal places precision
    return(GL)




# use pyvcf to read and parse the original VCF with the following formatting:
# GL's for SNPs that vary within this sample, GT:AD:GL e.g. 
# 0/1:7,2:28,0,167 means 0/1 is the called genotype; 7 ref allele & 2 alt allele counts; 28,0,167 are the genotype likelihoods
#vcf_reader = vcf.Reader(open("../data/maize282/vcf_AGPv4/subset_chr" + nChr + ".recode.vcf", 'r')) # small test file
vcf_reader = vcf.Reader(open('test8.recode.vcf', 'r'))
r = next(vcf_reader) # get first record (SNP)
GLs = []
for sample in r.samples:
    #print(sample['GL']) # for that record (SNP), prints genotype likelihoods for each sample
    GLs.extend(convertGL(sample['GL']))

print(GLs)


