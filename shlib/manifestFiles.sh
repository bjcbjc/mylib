

path=`pwd`
#manifestFiles.sh "file_pattern" destination_path

for fn in `ls $path/$1`
do
    echo -e `date +%Y"-"%m"-"%d`"\t"$fn"\t"$2"\t"bjchen"\t"manual
done

#2014-06-23      /data/analysis/SnyderA/Project_SNY_01395_RNA/Results/DESeq_normalized_count_matrix.txt  /data/delivery/SnyderA/data/Project_SNY_01395_RNA/Results       hmgeiger        manual
#2014-06-23      /data/analysis/SnyderA/Project_SNY_01395_RNA/Results/featureCounts_count_matrix.txt     /data/delivery/SnyderA/data/Project_SNY_01395_RNA/Results       hmgeiger        manual

#date +%Y"-"%m"-"%d | tr '\n' '\t' > /data/delivery_management/manifests/${person_name}_${project_name}_deliver_results.txt

#echo /data/analysis/${person_name}/${project_name}/Results/DESeq_normalized_count_matrix.txt | tr '\n' '\t' >> /data/delivery_management/manifests/${person_name}_${project_name}_deliver_results.txt

#echo "/data/delivery/${person_name}/data/${project_name}/Results hmgeiger manual" | awk '{ OFS="\t"}{$1=$1;print $0}' >> /data/delivery_management/manifests/${person_name}_${project_name}_deliver_results.txt




