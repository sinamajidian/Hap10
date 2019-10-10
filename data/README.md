

## A test for Fast mode.


You can test the pipeline using the provided data in this folder. We suppose that you have a BAM and VCF file. Also, you've built the extract_poly and SDhaP package .



```

cd data

k=3
mean_10x_molecule_length=50   


../utilities/break_vcf.sh  var.vcf var_break.vcf

cat var_break.vcf  | grep -v "1/1/1" | grep -v "0/0/0" grep -v "/2" > var_het.vcf  

./extract_poly/build/extractHAIRS --10X 1 --bam out_longranger/outs/possorted_bam.bam --VCF var_het.vcf --out unlinked_fragment_file

grep -v "NULL"  unlinked_fragment_file > unlinked_fragment_file_filtered  # removing those reads without barcode

python3 ../utilities/LinkFragments_brcd_based.py  unlinked_fragment_file_filtered frag.txt


python3 ../utilities/splitter.py frag.txt $mean_10x_molecule_length frag_sp.txt 

python3 ../utilities/extract_scc.py frag_sp.txt scc ./out



for i in {0..1}; do

	python2 ../utilities/FragmentPoly.py -f frag_sp0_${i}.txt  -o frag0_${i}_sd.txt -x SDhaP  
	hap_poly frag0_${i}_sd.txt out.hap $k 

done

```




