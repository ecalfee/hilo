#!/bin/bash
cd ~/hilo/data/landraces_fromLi/original

for i in $(awk '{print $1}' ../alloMaizeInclude.list); \
do echo $i; \
iget -K /iplant/home/lilepisorus/landraces_AGPv4/$i.bam; \
done
