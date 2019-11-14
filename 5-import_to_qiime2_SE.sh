working_dir=/home/ngunjiri1/Working_dir/data/June2019/sequences/demux
raw_data=/home/ngunjiri1/Working_dir/data/June2019/sequences/Primer-Filtered/Merged/filtered_by_length

rm -rf $working_dir
mkdir $working_dir

cd $raw_data

gzip *

qiime tools import \
  --type 'SampleData[SequencesWithQuality]'  \
  --input-path $raw_data \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path $working_dir/merged-sequences.qza