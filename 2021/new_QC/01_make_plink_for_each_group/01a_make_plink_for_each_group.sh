
cd /home/common/projects/ovine_selection/New_QC_sheep_2021/data/plink_files

plink --bfile /home/common/projects/ovine_selection/Data/plink_files/raw_plink_192_TOP/merged_top_v3 --keep /home/common/projects/ovine_selection/New_QC_sheep_2021/data/argali_id.txt --autosome-num 26 --make-bed --out argali
plink --bfile /home/common/projects/ovine_selection/Data/plink_files/raw_plink_192_TOP/merged_top_v3 --keep /home/common/projects/ovine_selection/New_QC_sheep_2021/data/romanovka_id.txt --autosome-num 26 --make-bed --out romanovka
plink --bfile /home/common/projects/ovine_selection/Data/plink_files/raw_plink_192_TOP/merged_top_v3 --keep /home/common/projects/ovine_selection/New_QC_sheep_2021/data/katahdin_id.txt --autosome-num 26 --make-bed --out katahdin
plink --bfile /home/common/projects/ovine_selection/Data/plink_files/raw_plink_192_TOP/merged_top_v3 --keep /home/common/projects/ovine_selection/New_QC_sheep_2021/data/hybrids_id.txt --autosome-num 26 --make-bed --out hybrids
