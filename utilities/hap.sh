#!/bin/bash



k=3	# ploidy level
m=50    # mean_10x_molecule_length
bam=$1
vcf=$2

./extract_poly/build/extractHAIRS --10X 1 --bam $bam --VCF $vcf  --out unlinked_fragments

python3 utilities/LinkFragments_brcd_based.py  unlinked_fragments frag.txt


python3 utilities/splitter.py frag.txt $m frag_sp.txt

python3 utilities/extract_scc.py frag_sp.txt scc ./out

ls | grep "frag_scc" | head -n -1 > list_frags_scc.txt



while read frag_i; do 

	python2 utilities/FragmentPoly.py -f ${frag_i}  -o ${frag_i}_sd.txt -x SDhaP
	hap_poly ${frag_i}_sd.txt  ${frag_i}.hap $k

	cat ${frag_i}.hap >> haplotype.hap
	
done < list_frags_scc.txt



python2 utilities/ConvertAllelesSDhaP.py -p haplotype.hap -o haplotype_with_genomic_position.hap -v $vcf  
python utilities/hap2vcf.py haplotype_with_genomic_position.hap $vcf

