Hufford 2012 Domestication Genes downloaded from Supplementary table 6:

 Letter | Published: 03 June 2012
Comparative population genomics of maize domestication and improvement

    Matthew B Hufford, Xun Xu, Joost van Heerwaarden, Tanja Pyhäjärvi, Jer-Ming Chia, Reed A Cartwright, Robert J Elshire, Jeffrey C Glaubitz, Kate E Guill, Shawn M Kaeppler, Jinsheng Lai, Peter L Morrell, Laura M Shannon, Chi Song, Nathan M Springer, Ruth A Swanson-Wagner, Peter Tiffin, Jun Wang, Gengyun Zhang, John Doebley, Michael D McMullen, Doreen Ware, Edward S Buckler, Shuang Yang & Jeffrey Ross-Ibarra

Nature Genetics volume 44, pages 808–811 (2012)

This paper used a comparison of XP-CLR cross-population composite likelihood ratio statistic (Chen 2010) 
to identify regions of the genome with elevated divergence between parviglumis and maize landraces
extending over large regions, consistent with fast rises in frequency (ie selection) during domestication
Candidate genes were identified within these elevated regions. One concern is that adaptive introgression
or even neutral introgression, from mexicana might also create elevated XP-CLR (not sure).

I have translated the gene names from Hufford's original paper to the names used for the maize v4 genome assembly:
https://www.maizegdb.org/gene_center/gene#
copy-and-paste into "Translate Gene Model IDs" tool. Selected 'translate to' Zm00001d.2.
Saved output file as gene_model_translation_to_APGv4.txt
Saved just list of gene names on APGv4 to 
domestication$ cut -f 3 gene_model_translation_to_APGv4.txt | uniq | awk '$0 != "" && $0 != "Zm00001d.2" {print $0}' | \
head > genes_on_APGv4_Zm00001d.2_gene_model.list
# I think the third column is the one that can be found in APGv4 gff3 file.
# My goal is to make 2 bed files: all the domestication genes and all non-domestication genes
# first I limit the gff3 to just genes (~45K):
data/refMaize/geneAnnotations$ cat Zea_mays.B73_RefGen_v4.41.chr.gff3 | awk -F" " '$3 ~ /gene/ {print $0}' > Zea_mays.B73_RefGen_v4.41.chr.genes.only.gff3
# then I get all genes with a match in Zm00001d.2 gene set (which I believe should all be in APG4 annotations)
data/domestication$ awk '$3 != "" {print $3}' gene_model_translation_to_APGv4.txt | tail -n +2 > gene_model_translation_to_APGv4_Zm00001d.2_hits.txt
# now I can look at ancestry calls for these genes and for other genes in R. But I may have to create windows around the genes for ancestry
# because genes are so short it's unlikely they will have an ancestry call within...but I could average the 2 closest ancestry calls on either side or another window approach
# first I want to get ancestry at all genes and then use permutations to see if domestication genes look different. 
# I could also look at depth of coverage across domestication genes.
