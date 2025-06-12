cd ~/rds/rds-rb643-ukbiobank2/Data_Users/yh464/params/ref/1000g_by_eth
for eth in eas amr afr eur sas; do
  for chrom in {1..22}; do 
    plink --bfile $eth/all_chrs --chr $chrom --make-bed --out $eth/chr$chrom && ./lava-ldblk/ldblock $eth/chr$chrom -out $eth/chr$chrom & 
  done
done
