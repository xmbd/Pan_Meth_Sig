#!~/bin/bash
######## data prepared
#### Create a 100bp area around the DMPs of each signature
for i in *txt;do awk '{print $1"\t"$2-99"\t"$3+99"\t"$4"\t"$5"\t"$6}' $i|sort -k1,3V -k2|uniq|sed 's/chr//g' >> $i.bed;done
#### create a backgroud file
cat ~/100bp/backgroud_file_450k.txt|awk '{print $1"\t"$2-99"\t"$3+99"\t"$4"\t"$5"\t"$6}'|sort -k1,3V -k2|uniq|sed 's/chr//g' > ~/100bp//methy_450k_bg.bed
#### homer analysis
for i in *txt.bed.bed;do findMotifsGenome.pl $i ~/reference/Genome/hg19/hg19.fa result/$i -bg ~/100bp/methy_450k_bg.bed;done


