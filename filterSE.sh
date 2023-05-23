while read -r line;
do
echo ${line}
awk 'BEGIN{IFS=FS="\t";while(getline < "rmats_SE.txt"){d[$1]=$1}}NR==1{i=1;name[i]=$i;for(i=2;i<=NF;i++){split($i,n,"-");id=substr(n[4],1,2);if(id>=11){name[i]=$i}}}$1 in d{out=$1;for(i=2;i<=NF;i=i+1){if(i in name){out=out"\t"$i}}print out}' rmats_psi/tcga_${line}_psi_rMATS_v32_filter.tsv | sort | uniq > filtered_psi/${line}_psi_normal.txt
awk 'BEGIN{IFS=FS="\t";while(getline < "rmats_SE.txt"){d[$1]=$1}}NR==1{i=1;name[i]=$i;for(i=2;i<=NF;i++){split($i,n,"-");id=substr(n[4],1,2);if(id<11){name[i]=$i}}}$1 in d{out=$1;for(i=2;i<=NF;i=i+1){if(i in name){out=out"\t"$i}}print out}' rmats_psi/tcga_${line}_psi_rMATS_v32_filter.tsv | sort | uniq > filtered_psi/${line}_psi.txt
done < sel_TCGA_normal.list 

while IFS=$'\t' read -r -a line;
do
echo ${line[0]}
sbatch -o ${line[0]}_pos.out -e ${line[0]}_pos.err -p hpc -N 1 -n 1 -c 1 --wrap="Rscript get_sig_eve_tcga.R ${line[0]}"
done < sel_TCGA_normal.list


while read -r line;
do
echo ${line}
grep "chr10@+@104008356^104018160-@104008356^104010816-104010908^104018160-" filtered_psi/${line}_psi.txt > SLK_psi/${line}_psi.txt
grep "chr10@+@104008356^104018160-@104008356^104010816-104010908^104018160-" filtered_psi/${line}_psi_normal.txt > SLK_psi/${line}_psi_normal.txt
done < sel_TCGA_normal.list 

while read -r line;
do
echo ${line}
awk 'BEGIN{IFS=FS="\t"}NR==1{i=1;name[i]=$i;for(i=2;i<=NF;i++){split($i,n,"-");id=substr(n[4],1,2);if(id<11){name[i]=$i}}}$1~"23543"{out=$1;for(i=2;i<=NF;i=i+1){if(i in name){out=out"\t"$i}}print out}' gene_counts/${line}_gene_quantification.txt | sort | uniq > RBFOX2_exp/${line}_count.txt
awk 'BEGIN{IFS=FS="\t"}NR==1{i=1;name[i]=$i;for(i=2;i<=NF;i++){split($i,n,"-");id=substr(n[4],1,2);if(id>=11){name[i]=$i}}}$1~"23543"{out=$1;for(i=2;i<=NF;i=i+1){if(i in name){out=out"\t"$i}}print out}' gene_counts/${line}_gene_quantification.txt | sort | uniq > RBFOX2_exp/${line}_count_normal.txt
done < sel_TCGA_normal.list
