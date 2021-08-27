#!/bin/bash



k=3	# ploidy level
m=50    # mean_10x_molecule_length
bam=$1
vcf=$2



extractHAIRS --10X 1 --mbq 4  --mmq 4  --bam $bam --VCF $vcf --out unlinked_fragment_file

grep -v "NULL" unlinked_fragment_file > unlinked_fragment_file_filtered 


echo "Raw fragment file is generated."

python3 ../utilities/LinkFragments_brcd_based.py unlinked_fragment_file_filtered frag_linked.txt

echo "Molecule-specific fragment file is generated."

python3 ../utilities/splitter.py frag_linked.txt  $vcf $m


echo "Barcode-specific fragment file is generated."

python3 ../utilities/extract_scc.py frag_linked_sp.txt scc ./frag_scc

echo "Several fragment files including strongly connected components are generated."


ls | grep "frag_scc" | head -n -1 > list_frags_scc.txt

while read frag_i; do 

	python2 ../utilities/FragmentPoly.py -f ${frag_i} -o ${frag_i}_sd -x SDhaP  
	hap_poly ${frag_i}_sd ${frag_i}.hap $k
	cat ${frag_i}.hap >> haplotype.hap
	
done < list_frags_scc.txt

python2 ../utilities/ConvertAllelesSDhaP.py -p haplotype.hap -o haplotype_with_genomic_position.hap -v  $vcf  


python2 ../utilities/hap2vcf.py haplotype_with_genomic_position.hap $vcf
