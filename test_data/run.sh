#!/bin/bash



k=3	# ploidy level
m=50    # mean_10x_molecule_length
bam=$1
vcf=$2


../extract_poly/build/extractHAIRS --10X 1 --mbq 4  --mmq 4  --bam $bam --VCF $vcf --out unlinked_fragment_file

grep -v "NULL" unlinked_fragment_file > unlinked_fragment_file_filtered 

python3 ../utilities/LinkFragments_brcd_based.py unlinked_fragment_file_filtered frag.txt

python3 ../utilities/splitter.py frag.txt  $vcf $m

python3 ../utilities/extract_scc.py frag_sp.txt scc ./frag_scc


for i in {0..2}; do

	python2 ../utilities/FragmentPoly.py -f frag_scc0_${i}.txt  -o frag_scc0_${i}_sd.txt -x SDhaP  
	../utilities/hap_poly frag_scc0_${i}_sd.txt scc0_${i}.hap $k
	cat scc0_${i}.hap >> haplotype.hap
	
done

python2 ../utilities/ConvertAllelesSDhaP.py -p haplotype.hap -o haplotype_with_genomic_position.hap -v  $vcf  


python2 ../utilities/hap2vcf.py haplotype_with_genomic_position.hap $vcf
