#!/bin/bash
#example: sh vcf_comparison.sh sample_name (NIDDL_1234)
identifier=$1
vdir=vcf_comparison_dir

bgzip "$identifier"_snpeff_lofreq.vcf
bgzip "$identifier"_snpeff_allvariants.vcf
bcftools index "$identifier"_snpeff_lofreq.vcf.gz
bcftools index "$identifier"_snpeff_allvariants.vcf.gz
bcftools isec -c none -p ${vdir} "$identifier"_snpeff_lofreq.vcf.gz "$identifier"_snpeff_allvariants.vcf.gz 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AF\t%QUAL\t%DP4\t%ANN\n' ${vdir}/0002.vcf > ${vdir}/"$identifier"_comparison_lofreq.tsv
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%REF_DP][\t%REF_QUAL][\t%ALT_DP][\t%ALT_QUAL][\t%ALT_FREQ][\t%ANN][\t%SAMPLE=%GT]\n' ${vdir}/0003.vcf > ${vdir}/"$identifier"_comparison_iVar.tsv
awk -v sample="$identifier" -i inplace -F'\t' -vOFS='\t' ' { gsub("," , "\t" , $8); gsub(",","\t",$9);gsub(",","\t",$10);gsub(",","\t",$11); print sample "\t" $0}' ${vdir}/"$identifier"_comparison_lofreq.tsv 
awk -v sample="$identifier" -i inplace '{print sample "\t" $0}' ${vdir}/"$identifier"_comparison_iVar.tsv 
echo -e "SAMPLE\tREF_GENOME\tPOS\tREF\tALT\tDP\tAF\tQUAL\tDP_REF_FW\tDP_REF_RV\tDP_ALT_FW\tDP_ALT_RV\tANN\tREF_GENOME\tPOS\tREF\tALT\tREF_DP\tREF_QUAL\tALT_DP\tALT_QUAL\tALT_FREQ\tANN\tSAMPLE" > ${vdir}/"$identifier"_comparison_combined.tsv 
paste -d' ' ${vdir}/"$identifier"_comparison_lofreq.tsv ${vdir}/"$identifier"_comparison_iVar.tsv >> ${vdir}/"$identifier"_comparison_combined.tsv
