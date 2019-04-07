program: fqtool  
version: 0.0.0  
updated: 18:51:01 Apr  7 2019  
Usage: fqtool [OPTIONS]  

|Options                                                                        |Explanations
|-------------------------------------------------------------------------------|---------------------------------
|  -h,--help                                                                    | Print this help message and exit
|  -i,--in1 FILE REQUIRED                                                       | read1 input file name
|  -o,--out1 TEXT REQUIRED                                                      | read1 output file name
|  -I,--in2 FILE Needs: --in1 Excludes: --inverleaved_in                        | read2 input file name
|  -O,--out2 TEXT Needs: --in2                                                  | read2 output file name
|  -m,--merge                                                                   | merge overlapped pair of read into one
|  -d,--discard_unmerged Needs: --merge                                         | discard unmerged reads
|  --phred64                                                                    | input fastq quality is in phred64 mode
|  -z,--compress_level INT in [1 - 9]                                           | compression level for gzip output
|  --inverleaved_in Excludes: --in2                                             | input is an interleaved FASTQ
|  --notoverwrite                                                               | do not overwriting existing output
|  -v,--verbose                                                                 | output verbose log information
|  -A,--enable_adapter_trimming                                                 | enable adapter trimming
|  --adapter_seqr1 TEXT Needs: --enable_adapter_trimming                        | adapter for read1
|  --adapter_seqr2 TEXT Needs: --enable_adapter_trimming                        | adapter for read2
|  --detect_pe_adapter Needs: --in2                                             | enable detect adapter for PE reads
|  -f,--trim_front1 INT in [0 - 1000]=0                                         | #bases trimmed in front for read1
|  -t,--trim_tail1 INT in [0 - 1000]=0                                          | #bases trimmed in tail for read1
|  -b,--max_len1 INT in [0 - 1000]=0                                            | maximum length keep for read1 after trim tail
|  -F,--trim_front2 INT in [0 - 1000]=0                                         | #bases trimmed in front for read2
|  -T,--trim_tail2 INT in [0 - 1000]=0                                          | #bases trimmed in tail for read2
|  -B,--max_len2 INT in [0 - 1000]=0                                            | maximum length keep for read2 after trim tail
|  -g,--trim_polyg                                                              | enable polyG trimming
|  --polyg_min_len INT=10 Needs: --trim_polyg                                   | minimum length to detect polyG in the read tail
|  -x,--trim_polyx                                                              | enable polyX trimming
|  --polyx_min_len INT=10 Needs: --trim_polyx                                   | minimum length to detect polyG in the read tail
|  --cut_front                                                                  | slide window from 5' to 3', drop the bases in the window if its mean quality < threshold, stop otherwise.
|  --cut_tail                                                                   | slide window from 3' to 5', drop the bases in the window if its mean quality < threshold, stop otherwise.
|  --cut_right                                                                  | slide window from 5' to 3', if a window with mean quality < threshold found, drop bases in it and the right part then stop.
|  --cut_window_size INT in [0 - 1000]=4                                        | the window size option shared by cut_front, cut_tail or cut_sliding
|  --cut_mean_quality INT in [1 - 36]=20                                        | the mean quality requirement option shared by cut_front, cut_tail or cut_sliding.
|  --cut_front_window_size INT in [0 - 1000]=4 Needs: --cut_front               | the window size option of cut_front
|  --cut_tail_window_size INT in [0 - 1000]=4 Needs: --cut_tail                 | the window size option of cut_tail
|  --cut_right_window_size INT in [0 - 1000]=4 Needs: --cut_right               | the window size option of cut_right
|  --cut_front_mean_quality INT in [1 - 36]=20 Needs: --cut_front               | the mean quality option of cut_front
|  --cut_tail_mean_quality INT in [1 - 36]=20 Needs: --cut_tail                 | the mean quality option of cut_tail
|  --cut_right_mean_quality INT in [1 - 36]=20 Needs: --cut_tail                | the mean quality option of cut_right
|  --enable_quality_filtering                                                   | enable quality filtering
|  --qualified_quality_phred INT=0 Needs: --enable_quality_filtering            | the minimum quality value that a base is qualified
|  --unqualified_base_limit INT=40 Needs: --enable_quality_filtering            | maximum low quality bases allowed in one read
|  --n_base_limit INT=5 Needs: --enable_quality_filtering                       | maximum N bases allowed in one read
|  --enable_length_filter                                                       | enable length filter
|  --minimum_length INT in [0 - 1000]=15 Needs: --enable_length_filter          | minimum length required for a read
|  --maximum_length INT in [0 - 1000]=0 Needs: --enable_length_filter           | maximum length allowed for a read
|  --enabel_lowcomplexity_filter                                                | enable low complexity filter
|  --minimum_complexity INT in [0 - 1]=0.3 Needs: --enabel_lowcomplexity_filter | minimum complexity required for a read
|  --filter_by_index                                                            | enable index filtering
|  --filter_index1 FILE Needs: --filter_by_index                                | file contains a list of barcodes of index1 to be filtered out, one barcode per line
|  --filter_index2 FILE Needs: --filter_by_index                                | file contains a list of barcodes of index2 to be filtered out, one barcode per line
|  --filter_index_threshold INT in [0 - 10]=0 Needs: --filter_by_index          | allowed difference of index barcode for index filtering
|  --enable_base_correction                                                     | enable base correction in overlapped regions of PE reads
|  --overlap_len_required INT in [0 - 1000]=0                                   | minimum overlapped lenfth needed for overlap analysis based adapter trimming and correction
|  --overlap_diff_limit INT in [0 - 10]=5                                       | maximum difference allowed for overlapped region for overlap analysis based adapter trimming and correction
|  --enable_umi_processing                                                      | enable unique molecular identifier (UMI) preprocessing
|  --umi_loc INT in [1 - 6]=0 Needs: --enable_umi_processing                    | 0:none, 1:index1, 2:index2, 3:read1, 4:read2, 5:per_index, 6:per_read
|  --umi_len INT in [0 - 1000]=0 Needs: --enable_umi_processing                 | if the UMI is in read1/read2, its length should be provided
|  --umi_prefix TEXT Needs: --enable_umi_processing                             | if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI)
|  --umi_skip INT in [0 - 1000]=0 Needs: --enable_umi_processing                | if the UMI is in read1/read2, fastp can skip several bases following UMI
|  --enable_overrepana                                                          | enable overrepresented sequence analysis.
|  --overrepana_sampling INT in [1 - 10000]=20 Needs: --enable_overrepana       | one in --overrepana_sampling will be computed for overrepresentation analysis
|  --json TEXT=report.json                                                      | the json format report file name
|  --html TEXT=report.html                                                      | the html format report file name
|  --title TEXT=Fastq Report                                                    | report title
|  --thread INT in [1 - 16]=1                                                   | worker thread number
|  --split_by_file_number Excludes: --split_by_lines                            | split output by limiting total split file number
|  --file_number INT Needs: --split_by_file_number                              | total split output file number
|  --split_by_lines Excludes: --split_by_file_number                            | split output by limiting lines of each file
|  --file_lines UINT Needs: --split_by_lines                                    | split output file line limit
|  --split_prefix_digits INT in [1 - 10]=0                                      | the digits for sequential output

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
