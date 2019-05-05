#ifndef OPTIONS_H
#define OPTIONS_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "util.h"

/** struct to store buffer options */
struct BufferSizeOptions{
    size_t maxPacksInReadPackRepo; ///< max number of ReadPacks a ReadPackRepository can hold
    size_t maxReadsInPack;         ///< max number of reads a ReadPack can hold
    size_t maxPacksInMemory;       ///< max number of ReadPacks in memory allowed
    BufferSizeOptions(){
        maxPacksInReadPackRepo = 1000;
        maxReadsInPack = 1000;
        maxPacksInMemory = 500;
    }
};

/** struct to store merge options */
struct MergePEReadsOptions{
    bool enabled;           ///< enable merge overlapped pairend reads into one
    bool discardUnmerged;   ///< discard unmerged reads
    std::string out;        ///< output file name
    /** construct a MergePEReads object to store merge options */
    MergePEReadsOptions(){
        enabled = false;
        discardUnmerged = false;
    }
};

/** struct to store PolyG trimming options */
struct PolyGTrimmerOptions{
    bool enabled;       ///< enable PolyG trimming
    int minLen;         ///< minimum length needed to compare to get a proper polyG tail(total mismatch < 5 || each8basemismatches <= 1)
    /** construct a PolyGTrimmerOptions and set default values */
    PolyGTrimmerOptions(){
        enabled = false;
        minLen = 10;
    }
};

/** struct to store PolyX trimming options */
struct PolyXTrimmerOptions{
    bool enabled;       ///< enable PolyX trimming
    int minLen;         ///< minimum length needed to compare to get a proper polyX tail(total mismatch < 5 || each8basemismatches <= 1)
    /** construct a PolyXTrimmerOptions and set default values */
    PolyXTrimmerOptions(){    
        enabled = false;
        minLen = 10;
    }
};

/** struct to store umi process options */
struct UMIOptions{
    bool enabled;          ///< enable umi process if true
    int location;          ///< umi locations, refer to the constexpr defined
    int length;            ///< umi length
    int skip;              ///< number of bases to skip after umi cut from read
    bool dropOtherComment; ///< drop other comment string in header
    /** construct a UMIOptions object and set default values */
    UMIOptions(){
        enabled = false;
        location = 0;
        length = 0;
        skip = 0;
        dropOtherComment = false;
    }
};


/** struct to store duplication analysis options */
struct DuplicationAnalysisOptions{
    bool enabled;   ///< enable duplication analysis if true
    int keylen;     ///< key length of read, should be less than (std::numeric_limits<size_t>::digits  - 1) / 2
    int histSize;   ///< hist length to do statistics
    /** construct a DuplicationAnalysisOptions and set default values */
    DuplicationAnalysisOptions(){
        enabled = true;
        keylen = 12;
        histSize = 32;
    }
};

/** struct to store various quality threshold dependent fastq read cut/trim options */
struct QualityCutOptions{
    bool enableFront;      ///< if true, sliding window from front(5'end) until average quality > minFrontQual and cut from the last non-N base of this window
    bool enableTail;       ///< if true, sliding window from tail(3'end) until average quality > minTailQual and cut from the last non-N base of this window
    bool enableRright;     ///< if true, sliding window from front(5'end) until average quality < minRightQual and cut before the first base with quality < minRightQual
    int qualityShared;     ///< shared low quality threshold for cutting, make setting default values easy
    int windowSizeShared;  ///< shared window size for sliding, make setting default values easy
    int qualityFront;      ///< minimum average window quality required to stop sliding cut from front(5'end)
    int qualityTail;       ///< minimum average window quality required to stop sliding cut from tail(3'end)
    int qualityRight;      ///< minimum average window quality required to stop sliding detect a low quality window from front(5'end)
    int windowSizeFront;   ///< window size for sliding cut from front
    int windowSizeTail;    ///< window size for sliding cut from tail
    int windowSizeRight;   ///< window size for sliding detect low quality window from front
    
    /** construct a QualityCutOptions object and set default values */
    QualityCutOptions(){
        enableFront = false; 
        enableTail = false;  
        enableRright = false;
        qualityShared = 20;       
        windowSizeShared = 4;    
        qualityFront = qualityShared;     
        qualityTail = qualityShared;         
        qualityRight = qualityShared;        
        windowSizeFront = windowSizeShared;     
        windowSizeTail = windowSizeShared;      
        windowSizeRight = windowSizeShared;
    }
};

/** struct to store index filtering options */
struct IndexFilterOptions{
    bool enabled;                         ///< enable index filtering if set true
    int threshold;                        ///< maximum different bases allowed for an external index hit blacklist successfully 
    std::string index1File;               ///< index1 file path 
    std::string index2File;               ///< index2 file path
    std::vector<std::string> blacklist1;  ///< read1 index blacklist
    std::vector<std::string> blacklist2;  ///< read2 index blacklist

    /** construct a IndexFilterOptions object and set default values */
    IndexFilterOptions(){
        enabled = false;
        threshold = 0;
    }
};

/** struct to store over representation sequence analyisis options */
struct OverrepresentedSequenceAnalysisOptions{
    bool enabled;                                     ///< enable over representation sequence analyisis
    int sampling;                                     ///< sampling frequence for ORA
    std::map<std::string, size_t> overRepSeqCountR1;  ///< over represented sequences count of read1
    std::map<std::string, size_t> overRepSeqCountR2;  ///< over represented sequences count of read2
    /** construct a OverrepresentedSequenceAnalysisOptions object and set default values */
    OverrepresentedSequenceAnalysisOptions(){
        enabled = false;
        sampling = 20;
    }
};

/** struct to store base correction options */
struct  CorrectionOptions {
    bool enabled; ///< enable base correction if true
    /** construct a CorrectionOptions object and set default values */
    CorrectionOptions() {
        enabled = false;
    }
};

/** struct to store low complexity filter options */
struct LowComplexityFilterOptions {
    bool enabled;     ///< enable low complexity filter if true
    double threshold; ///< threshold to test sequence complexity, below which will be dropped 
    /** construct a LowComplexityFilterOptions object and set default values */
    LowComplexityFilterOptions() {
        enabled = false;
        threshold = 0.3;
    }
};

/** struct to store read length filter options */
struct ReadLengthFilterOptions{
    bool enabled;         ///< enable read length filter if true
    int minReadLength;    ///< if read_length < minReadLength, then this read will be discarded
    int maxReadLength;    ///< if read_length > maxReadLength, then this read will be discarded
    /** construct a ReadLengthFilterOptions object and set default values */
    ReadLengthFilterOptions(){
        enabled = false;
        minReadLength = 15;
        // 0 for no limination
        maxReadLength = 0; 
    }
};

/** struct to store quality filter options */
struct QualityFilterOptions{
    bool enabled;               ///< quality filter enabled if true
    char lowQualityLimit;       ///< if a base's quality < lowQualityLimit, hen it's considered as a low quality base
    int lowQualityBaseLimit;    ///< if low quality bases number > lowQualityBaseLimit, then discard this read
    int nBaseLimit;             ///< if N bases number > nBaseLimit, then discard this read
    double lowQualityRatio;     ///< if a read has bases with quality < lowQualityLimit more than this ratio, it will be dropped
    double averageQualityLimit; ///< if a read has averaeg quality less then averageQualityLimit, it will be dropped
    /** construct a QualityFilterOptions object and set default values */
    QualityFilterOptions(){
        enabled = true;
        // '0' == Q15
        lowQualityLimit = 48;
        lowQualityBaseLimit = 40;
        nBaseLimit = 5;
        lowQualityRatio = 0.30;
        averageQualityLimit = 0.0;
    }
};

/** struct to store adapter detect/trimming options */
struct AdapterOptions{
    bool cutable;                     ///< the adpter of reads can be cut if true
    bool enableTriming;               ///< enable index trimming if true
    bool enableDetectForPE;           ///< enable auto detection of index for pair end reads if true
    bool adapterSeqR1Provided;        ///< adapter sequence for read1 is provided externally if true
    bool adapterSeqR2Provided;        ///< adapter sequence for read2 is provided externally if true
    std::string inputAdapterSeqR1;    ///< adapter sequence for read1 provided externally
    std::string inputAdapterSeqR2;    ///< adapter sequence for read2 provided externally
    std::string detectedAdapterSeqR1; ///< adapter sequence for read1 auto detected 
    std::string detectedAdapterSeqR2; ///< adapter sequence for read2 auto detected
    double reportThreshold;           ///< adapter sequence trim count rate more than this value will be reported
    /** construct a AdapterOptions object and set default values */
    AdapterOptions(){
        cutable = false;
        enableTriming = true;
        adapterSeqR1Provided = false;
        adapterSeqR2Provided = false;
        enableDetectForPE = true;
        reportThreshold = 0.01;
    }
};

/** struct to store force trimming options */
struct ForceTrimOptions{
    int front1;  ///< first cycles trimmed for read1
    int tail1;   ///< last cycles trimmed for read1
    int front2;  ///< first cycles trimmed for read2
    int tail2;   ///< last cycles trimmed for read2
    int maxLen1; ///< maximum length of read1
    int maxLen2; ///< maximum length of read2
    /** construct a ForceTrimOptions object and set default values */
    ForceTrimOptions(){
        front1 = 0;
        tail1 = 0;
        front2 = 0;
        tail2 = 0;
        maxLen1 = 0;
        maxLen2 = 0;
    }
};

/** struct to store output file split options */
struct SplitOptions{
    bool enabled;        ///< enable output file split
    int number;          ///< split file numbers of a file
    size_t size;         ///< number of lines of each split file
    int digits;          ///< digits number of split filename prefix, e.g 0001 means 4 digits
    bool needEvaluation; ///< need evaluation if true
    bool byFileNumber;   ///< split by file number
    bool byFileLines;    ///< split by file lines
    /** construct a SplitOptions object and set default values */
    SplitOptions(){
        enabled = false;
        needEvaluation = false;
        number = 0;
        size = 0;
        digits = 4;
        byFileLines = false;
        byFileLines = false;
    }
};

/** struct to store Kmer analysis options */
struct KmerOptions{
    bool enabled;       ///< enable Kmer analysis if true
    int kmerLen;        ///< Kmer length to calculate 
    /** construct a KmerOptions object and set default values */
    KmerOptions(){
        enabled = false;
        kmerLen = 0;
    }
};

/** struct to store estimatation options */
struct EstimateOptions{
    int seqLen1;          ///< estimated read1 length
    int seqLen2;          ///< estimated read2 length
    int readsNum;         ///< estimated total read number
    bool twoColorSystem;  ///< estimated read from two color system if true
    std::string adapter;  ///< estimated adapter sequence
    bool illuminaAdapter; ///< estimated adapter sequnce is from illumina
    bool estimated;       ///< estimation work done already if true
    /** construct an EstimateOptions object and set default values */
    EstimateOptions(){
        seqLen1 = 151;
        seqLen2 = 151;
        readsNum = 0;
        twoColorSystem = false;
        adapter = "";
        illuminaAdapter = false;
        estimated = false;
    }
};

/** struct to hold various option structs and interface  options together */
struct Options{
    std::string version;          ///< program version
    std::string compile;          ///< program compile time
    std::string in1;              ///< input read1 filename
    std::string in2;              ///< input read2 filename
    std::string out1;             ///< output read1 filename
    std::string out2;             ///< output read2 filename
    std::string unpaired1;        ///< output unpaired read1 filename
    std::string unpaired2;        ///< output unpaired read2 filename
    std::string failedOut;        ///< output failed read1/2 filename
    std::string jsonFile;         ///< output json filename
    std::string htmlFile;         ///< output html report filename
    std::string reportTitle;      ///< html report title
    int digits;                   ///< number of digits for split filename prefix
    int compression;              ///< compression level for gz format output
    bool phred64;                 ///< the input file is using phred64 quality scoring if true 
    bool inputFromSTDIN;          ///< read from STDIN
    bool outputToSTDOUT;          ///< write to STDOUT
    bool interleavedInput;        ///< the input read1(in1) file is an interleaved PE fastq
    size_t readsToProcess;           ///< number of reads to process
    int thread;                   ///< number of threads to do paralel work
    int insertSizeMax;            ///< maximum value of insert size
    int overlapRequire;           ///< overlap region minimum length
    int overlapDiffLimit;         ///< overlap region maximum different bases allowed
    // submodule options
    ForceTrimOptions trim;                ///< ForceTrimOptions object
    QualityFilterOptions qualFilter;      ///< QualityFilterOptions object
    QualityCutOptions qualitycut;         ///< QualityCutOptions object 
    ReadLengthFilterOptions lengthFilter; ///< ReadLengthFilterOptions object 
    AdapterOptions adapter;               ///< AdapterOptions object 
    CorrectionOptions correction;         ///< CorrectionOptions object 
    OverrepresentedSequenceAnalysisOptions overRepAna; ///< OverrepresentedSequenceAnalysisOptions object
    LowComplexityFilterOptions complexityFilter;       ///< LowComplexityFilterOptions object
    IndexFilterOptions indexFilter;                    ///< IndexFilterOptions object
    SplitOptions split;                                ///< SplitOptions object
    KmerOptions kmer;                                  ///< KmerOptions object
    EstimateOptions est;                               ///< EstimateOptions object
    DuplicationAnalysisOptions duplicate;              ///< DuplicationAnalysisOptions object
    UMIOptions umi;                                    ///< UMIOptions object 
    PolyGTrimmerOptions polyGTrim;                     ///< PolyGTrimmerOptions object
    PolyXTrimmerOptions polyXTrim;                     ///< PolyXTrimmerOptions object
    MergePEReadsOptions mergePE;                       ///< MergePEReadsOptions object
    BufferSizeOptions bufSize;                         ///< BufferSizeOptions object
    std::mutex logmtx;                                 ///< mutex for logging
    
    /** Construct a Options object */
    Options();
    
    /** update options */
    void update();

    /** valid options */
    void validate();
    
    /** test input is paired or not
     * @return true if input is paired
     */
    bool isPaired();
    
    /** initialize index filter
     * @param blacklistFile1 index1 blacklist file
     * @param blacklistFile2 index2 blacklist file
     * @param threshold threshold for index and blacklist match
     */
    void initIndexFilter(const std::string& blacklistFile1, const std::string& blacklistFile2, int threshold = 0);
    
    /** get a vector of string from lines of a file
     * @param filename file path name
     * @return a vector of string from lines of a file
     */ 
    std::vector<std::string> makeListFromFileByLine(const std::string& filename);
};

#endif
