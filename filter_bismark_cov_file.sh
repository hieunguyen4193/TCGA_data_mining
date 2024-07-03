inputdir="/mnt/archiving/DATA_HIEUNHO/GW-bismark_cfDNA-highdepth_pair-end/06_methylation_extract";
files=$(ls ${inputdir}/*.cov | xargs -n 1 basename);
outputdir="/mnt/DATASM14/hieunho/hieu_project/metadata/ECD/outputdir";
ann450k_file="/datassd/hieunguyen/tmp/ann450kdf.txt"
for file in $files;do echo -e "Working on file " $file "\n" && bedtools intersect -a ${ann450k_file} -b ${inputdir}/${file} -wb > ${outputdir}/${cov_file%.cov*}.filtered450k.cov;done