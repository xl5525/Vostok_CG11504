~/gatk/gatk-4.6.1.0/gatk SelectVariants \
   -V brain5.new.b.vcf.gz \
   --exclude-non-variants \
   -O brain5.variant.vcf.gz


~/gatk/gatk-4.6.1.0/gatk SelectVariants \
   -V brain3.new.b.vcf.gz \
   --exclude-non-variants \
   -O brain3.variant.vcf.gz

~/gatk/gatk-4.6.1.0/gatk SelectVariants \
   -V cg11504_1.new.b.vcf.gz \
   --exclude-non-variants \
   -O cg11504_1.variant.vcf.gz

~/gatk/gatk-4.6.1.0/gatk SelectVariants \
   -V cg11504_2.new.b.vcf.gz \
   --exclude-non-variants \
   -O cg11504_2.variant.vcf.gz


bcftools isec -n=1 -w1 -p unique_cg11504_2 \
  cg11504_2.variant.vcf.gz \
  cg11504_1.variant.vcf.gz \
  brain3.variant.vcf.gz \
  brain5.variant.vcf.gz



~/gatk/gatk-4.6.1.0/gatk SelectVariants \
   -V cg11504_1.g.vcf \
   --exclude-non-variants \
   -O cg11504_1.g.variant.vcf

~/gatk/gatk-4.6.1.0/gatk SelectVariants \
   -V cg11504_2.g.vcf \
   --exclude-non-variants \
   -O cg11504_2.g.variant.vcf


~/gatk/gatk-4.6.1.0/gatk SelectVariants \
   -V brain3.g.vcf \
   --exclude-non-variants \
   -O brain3.g.variant.vcf

~/gatk/gatk-4.6.1.0/gatk SelectVariants \
   -V brain5.g.vcf \
   --exclude-non-variants \
   -O brain5.g.variant.vcf


bcftools isec -n=1 -p unique_g_out \
  brain3.g.variant.vcf.gz \
  cg11504_2.g.variant.vcf.gz \
  brain5.g.variant.vcf.gz \
  cg11504_1.g.variant.vcf.gz

awk 'BEGIN {OFS="\t"} {print $1, $2-1, $2, $3, $4, $5}' sites.txt > sites.bed





