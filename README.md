program: fqtool  
version: 0.0.0  
updated: 14:53:20 May 15 2019  
Usage: fqtool [OPTIONS]  

|Options                                                                       | Explanations
|------------------------------------------------------------------------------|---------------------------------
|  -h,--help                                                                   |  Print this help message and exit
|IO Options:  
|  -i FILE REQUIRED                                                            |  read1 input file name
|  -o TEXT REQUIRED                                                            |  read1 output file name
|  -I FILE Needs: -i Excludes: --in_fq_interleaved                             |  read2 input file name
|  -O TEXT Needs: -I                                                           |  read2 output file name
|  --unpaired_read1 TEXT                                                       |  output read1 whose mate failed QC
|  --unpaired_read2 TEXT                                                       |  output read2 whose mate failed QC
|  --failed_out TEXT                                                           |  output failed QC reads
|  --phred64                                                                   |  input fastq is phred64
|  -z INT in [1 - 9]                                                           |  gzip output compress level
|  --in_fq_interleaved Excludes: -I                                            |  input fastq interleaved
|Merge:  
|  -m Needs: -I Excludes: -s -S                                                |  merge overlapped readpair
|  --discard_unmerged Needs: -m                                                |  discard unmerged reads
|  --merge_output TEXT Needs: -m                                               |  merged output
|Duplication:
|  -d                                                                          |  enable duplication analysis
|  --dup_ana_key_len INT in [12 - 31]=12 Needs: -d                             |  duplication analysis key length
|  --dup_ana_hist_size INT in [1 - 10000]=32 Needs: -d                         |  duplicate analysis hist size
|Adapter:
|  -a                                                                          |  enable adapter trimming
|  --adapter_of_read1 TEXT Needs: -a                                           |  adapter of read1
|  --adapter_of_read2 TEXT Needs: -a                                           |  adapter of read2
|  --detect_pe_adapter Needs: -I                                               |  detect PE adapters
|Trim:
|  -f INT in [0 - 1000]=0                                                      |  bases trimmed in read1 front
|  -t INT in [0 - 1000]=0                                                      |  bases trimmed in read1 tail
|  -b INT in [0 - 1000]=0                                                      |  read1 max length allowed
|  -F INT in [0 - 1000]=0                                                      |  bases trimmed in read2 front
|  -T INT in [0 - 1000]=0                                                      |  #bases trimmed in read2 tail
|  -B INT in [0 - 1000]=0                                                      |  read2 max length allowed
|PolyX:
|  -g                                                                          |  enable polyG trim
|  --min_len_detect_polyG INT=10 Needs: -g                                     |  minimum length to detect polyG
|  -x                                                                          |  enable polyX trim
|  --min_len_detect_polyX INT=10 Needs: -x                                     |  minimum length to detect polyG
|Cut:
|  --enable_cut_front                                                          |  slide and drop from 5'->3'
|  --enable_cut_tail                                                           |  slide and drop from 3'->5'
|  --enable_cut_right                                                          |  slide from 5'->3' and drop window and right part
|  -W INT in [0 - 1000]=4                                                      |  window size for cut sliding
|  -M INT in [1 - 36]=20                                                       |  min mean quality to drop window/bases
|  --cut_front_window INT in [0 - 1000]=4 Needs: --enable_cut_front            |  window size to cut from 5''
|  --cut_tail_window INT in [0 - 1000] Needs: --enable_cut_tail                |  window size to cut from 3'
|  --cut_right_window INT in [0 - 1000]=4 Needs: --enable_cut_right            |  window size to cut right
|  --cut_front_mean_qual INT in [1 - 36]=20 Needs: --enable_cut_front          |  mean quality to cut from 5'
|  --cut_tail_mean_qual INT in [1 - 36] Needs: --enable_cut_tail               |  mean quality to cut from 3'
|  --cut_right_mean_qual INT in [1 - 36]=20 Needs: --enable_cut_tail           |  mean quality to cut right
|Qual:
|  -q                                                                          |  enable quality filter
|  -Q INT in [33 - 75]=0 Needs: -q                                             |  minimum ASCII Quality for qualified bases,
|  -U INT in [0 - 1]=0.3 Needs: -q                                             |  maximum low quality ratio allowed in one read
|  -N INT=5 Needs: -q                                                          |  maximum N bases allowed in one read
|  -e FLOAT Needs: -q                                                          |  average quality needed for one read
|Length:
|  -l                                                                          |  enable length filter
|  --min_length INT in [0 - 1000]=15 Needs: -l                                 |  min length required for a read
|  --max_length INT in [0 - 1000]=0 Needs: -l                                  |  max length allowed for a read
|Complexity:
|  -y                                                                          |  enable low complexity filter
|  -Y INT in [0 - 1]=0.3 Needs: -y                                             |  min complexity required for a read
|Index:
|  --enable_index_filter                                                       |  enable index filtering
|  --index1_file FILE Needs: --enable_index_filter                             |  index1 file to filter
|  --index2_file FILE Needs: --enable_index_filter                             |  index2 file to filetr
|  --max_diff_for_match INT in [0 - 10]=0 Needs: --enable_index_filter         |  max ed to validate index matcha
|Correction:
|  -c                                                                          |  enable base correction in PE reads
|  --min_overlap_len INT in [0 - 1000]=30                                      |  min overlap length needed for overlap analysis
|  --max_diff_for_overlap INT in [0 - 10]=5                                    |  max ed to validate overlap
|UMI:
|  -u                                                                          |  enable UMI preprocess
|  --umi_location INT in [1 - 6]=0 Needs: -u                                   |  0[none]1[index1]2[index2]3[read1]4[read2]5[perindex]6[perread]
|  --umi_length INT in [0 - 1000]=0 Needs: -u                                  |  umi length
|  --umi_skip_length INT in [0 - 1000]=0 Needs: -u                             |  bases to skip after umi
|  --umi_drop_comment Needs: -u                                                |  drop other comment information
|  --umi_not_trim Needs: -u                                                    |  do not trim reads
|ORA:
|  --ora                                                                       |  enable ORA
|  --ora_sample INT in [1 - 10000]=20 Needs: --ora                             |  ORA sampling steps
|KMer:
|  --kmer                                                                      |  enable kmer analysis
|  --kmer_length INT in [4 - 16]=0 Needs: --kmer                               |  kmer length to analysis
|Report:
|  -J TEXT=report.json                                                         |  json format report file
|  -H TEXT=report.html                                                         |  html format report file
|System:
|  -w INT in [1 - 16]=4                                                        |  worker thread number
|  --max_packs_in_repo INT in [1 - 1000000]=1000                               |  max packs in repo
|  --max_item_in_pack INT in [1 - 1000000]=100000                              |  max read/pairs in pack
|  --max_packs_in_mem INT in [1 - 1000000]=5                                   |  max packs in memory
|Split:
|  -s Excludes: -m -S                                                          |  split output by file number
|  --split_file_number INT Needs: -s                                           |  total split output file number
|  -S Excludes: -m -s                                                          |  max line of each output file
|  --splie_file_line UINT Needs: -S                                            |  split output file line limit
|  --digits_file_name INT in [1 - 10]=0                                        |  digits for sequential output filename

Installation  

1. clone repo  
`git clone https://github.com/vanNul/fqtool`  

2. compile  
`cd pipe`  
`./autogen.sh`  
`./configure --prefix=/path/to/install/dir/`  
`make`   
`make install`  

3. execute  
`/path/to/install/dir/bin/fqtool`  

4. test   
`fqtool -i ./testdata/r1.fq.gz -I ./testdata/r2.fq.gz -o r1.out.fq.gz -O r2.out.fq.gz-q --kmer --kmer_length 6 -d -a --detect_pe_adapter > run.log 2>&1` 

PS  
fqtool is modified from [fastp](https://github.com/OpenGene/fastp) with many enchancements, mainly listed below 
1. all ringbuffer are fixed and are real ring now, so fqtool will not crash with too many reads in input fastq
2. use json.hpp and ctml.hpp to write json and html respectly for further customized development with ease 
3. kmer analysis is disabled by default, and the kmer length used for analysis can be set 
4. buffer size arguments are externalized as well  
5. umi use standardized OX:Z tag and umi quality use standarzied BZ:Z to append as comments after readnames  
6. fastq readnumber estimation and overlap analysis fixedup  
