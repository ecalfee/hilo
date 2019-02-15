# These bams were obtained from Li Wang.
# Samples are individual landraces from across the Americas, analyzed in Wang 2017 https://doi.org/10.1186/s13059-017-1346-4
# The original analysis was done on v3 of the maize genome, but Li remapped all of them to v4 (official release, all contigs)
# One sample the remapped bam was corrupted and so I obtained the original v3 bam from Li:
# RIMMA0625_Andean.IndelRealigned.bam
# And will process v3 bam -> fastq -> map to v4 bam

# Note that RIMA0625 and RIMMA0438 had population labels mistakenly switched in the original genetics analysis  
# http://www.genetics.org/content/genetics/suppl/2015/06/15/genetics.115.178327.DC1/genetics.115.178327-4.pdf
# (see also mislabel here: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP065483)

# Li identified this error. RIMA0625 is from the Andes and RIMMA043 is from the Mexican highlands. This correction is reflected in the ID by population .txt files.
