

## A simple instruction for testing fast mode.


You can test the pipeline using the provided data (BAM and VCF files) in this folder.  we supposed that you've built the extract_poly and SDhaP package.   



```
chmod 755 run.sh
./run.sh possorted_bam.bam var_het.vcf

```

The final haplotype file is in both VCF format or a text `haplotype_with_genomic_position.hap`.


