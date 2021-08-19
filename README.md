# HILO project
Scripts related to the following publication:

Calfee E, Gates D, Lorant A, Perkins MT, Coop G, Ross-Ibarra J (2021). Selective sorting of ancestral introgression in maize and teosinte along an elevational cline. bioRxiv. https://www.biorxiv.org/content/10.1101/2021.03.05.434040v1

Analyses can be reproduced using snakemake (see Snakefile). Because the full pipeline is very long, scripts are divided into different sub pipelines, each with their own snakefile.

Where to find scripts for the main figures:
- Fig 1: map/
- Fig 2: global_ancestry/
- Fig 3: diversity/
- Fig 4: ancestry_by_r/
- Figs 5-7: ZAnc/

It is not recommended to run the full pipeline from scratch all at once. To use the main snakefile to reproduce results piece by piece, you can comment/uncomment subpipelines and outputs as you go along. A log of all commands run (including failed runs etc.) is in snake_commands.txt and, for earlier commands, older_snake_commands.txt.
