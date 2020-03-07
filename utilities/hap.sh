#!/bin/bash



k=3	# ploidy level
m=50    # mean_10x_molecule_length
bam=$1
vcf=$2

./extract_poly/build/extractHAIRS --10X 1 --bam $bam --VCF $vcf  --out unlinked_fragments

python3 utilities/LinkFragments_brcd_based.py  unlinked_fragments frag.txt


python3 utilities/splitter.py frag.txt $m frag_sp.txt

python3 utilities/extract_scc.py frag_sp.txt scc ./out



python2 utilities/FragmentPoly.py -f frag_sp.txt  -o frag_sd.txt -x SDhaP
./sdhap/hap_poly frag_sd.txt  out.hap 3
python2 utilities/ConvertAllelesSDhaP.py -p out.hap -o out_with_genomic_position.hap -v var_het.vcf  
python utilities/hap2vcf.py haplotype_with_genomic_position.hap $vcf

