Results of brief literature search for major known inversions:
These 4 inversions and their positions are based on LD patterns in 21 wild teosinte populations from Fig 3 in Pyhajarvi 2013:
Listed by them are other papers, sometimes with distinct methods, that also identified these inversions (see text near Fig. 3 for more information)
Also anything in old coordinates was converted to APGv4 genome assembly with this tool: http://ensembl.gramene.org/Zea_mays/Tools/AssemblyConverter?db=core
Inv1n 1:64389681-116086290  [Fang et al. 2012; Hufford et al. 2012] - APGv2 coordinates
Inv4m 4:168832447-182596678 [Mano et al. 2012] - APGv2 coordinates
Inv9e_or_d 9:0-29403252 [Ting 1964] -- citation & coordinates could be switched; I'm not sure which is 9e vs. 9d - APGv2 coordinates
Inv9d_or_e 9:118862807-143163957 - APGv2 coordinates

Separately from Pyhajarvi, in a large maize mapping population, Romero Navarro 2017 find
"a 6-Mb region on chromosome 3 beginning at 79 Mb. The 6-Mb region on chromosome 3 has similar segregation to that of Inv4,
and its increased LD suggests that it might be an inversion. In NAM populations, this putative inversion
and the centromere comprise a single quantitative trait locus (QTL) for flowering time" (Romero 2017)
Inv3flower 3:79000000-85000000 - GBS build 2.7 using TASSEL (this GBS build is aligned to APGv2 coordinates)

Additionally, the US-Chinese nested association mapping panel (NAM), Rodgers-Melnick 2015 PNAS study,
found 3 putative inversions with no crossovers in some family crosses (suggestive of segregating for an inversion),
one of which was consistent with a fixed inversion between B73 and the outgroup, sorghum
Inv1NAM 1:217900000-245500000 "Figure S13. The inversion of the B73 chromosome 1 segment 217.9-245.48 Mb relative to sorghum chromosome 1. In B73 x CML333 this region contains no
crossovers." (from mapping it goes up to 245.5) - APGv2 coordinates
Inv3NAM 3:167200000–176500000 (no crossovers in some families -- no found inversion with sorghum synteny) - APGv2 coordinates
Inv5NAM 5:177800000–194100000 (also no crossovers in some families -- no found inversion with sorghum synteny) - APGv2 coordinates
# I've put these coordinates in a file refMaize/inversions/knownInv.txt (all in v2 coordinates)

# All v2 coordinates were converted to v4 using gramene Zea Mays web tool Assembly Converter (http://ensembl.gramene.org/Zea_mays/Tools/AssemblyConverter). Then I used bedtools to merge all ranges closer than 1Mb
and identified the true inversion coordinates as the only merged range > 5Mb:
inversions/v4_coord$ for bed in $(ls *.bed); do bedtools sort -i $bed | bedtools merge -d 1000000 -i - | awk '{print $0"\t"$3-$2}' | awk -v b=${bed%.*} '$4>=5000000 {print b"\t"$0}'; done > ../knownInv_v4_coord.txt
