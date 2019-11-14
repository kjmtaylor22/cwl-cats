working_dir=/home/ngunjiri1/Working_dir/data/June2019/sequences/Primer-Filtered

raw_data=/home/ngunjiri1/Working_dir/data/June2019/sequences/Adapter_Trimmed

tag5='GAGTGCCAGCMGCCGCGGTAA'
tag5_R='ACGGACTACHVGGGTWTCTAAT'
      

rm -rf $working_dir
mkdir $working_dir
mkdir $working_dir/fastq

cd $raw_data


ls *.fastq > $raw_data/list.txt

while read R1
do read R2

echo "Working on $R1 $R2"

A="$( echo $R1 | cut -d'_' -f1 | sed s'/TOS//')"
B="$( echo $R1 | cut -d'_' -f2 | sed s'/S//')"
out_name_R1=$A"_"$B"_L001_R1_001"
out_name_R2=$A"_"$B"_L001_R2_001"
echo $out_name_R1

cutadapt -g $tag5 -G $tag5_R -o $working_dir/$out_name_R1".fastq"  -p $working_dir/$out_name_R2".fastq" $R1  $R2 -j 35 -e 0.2 -m 50


gzip $working_dir/$out_name_R1".fastq" 
gzip $working_dir/$out_name_R2".fastq"

mv $working_dir/*.gz $working_dir/fastq/
done< $raw_data/list.txt

