#ifndef FILTER_RESULT_H
#define FILTER_RESULT_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <json.hpp>
#include "common.h"
#include "options.h"
#include "stats.h"
#include "htmlutil.h"

/** Class to store/calculate read filter results */
class FilterResult{
    public:
        Options *mOptions;                                         ///< pointer to Options object
        bool mPaired;                                              ///< pairend filter results if true, else single end results
        size_t mCorrectedReads;                                    ///< number of reads got bases corrected
        size_t mCorrectedBases;                                    ///< number of bases got corrected
        size_t mMergedPairs;                                       ///< number of read pairs merged
        size_t mFilterReadStats[COMMONCONST::FILTER_RESULT_TYPES]; ///< reads filter status accumulation array
        size_t mTrimmedAdapterReads;                               ///< number of reads got adapter trimmed
        size_t mTrimmedAdapterBases;                               ///< number of bases got trimmed when do adapter trimming
        size_t* mTrimmedPolyXReads;                                ///< number of reads got trimmed when do polyx trimming
        size_t* mTrimmedPolyXBases;                                ///< number of bases got trimmed when do polyx trimming
        std::map<std::string, size_t> mAdapter1Count;              ///< read1 adapter trimmed count map
        std::map<std::string, size_t> mAdapter2Count;              ///< read2 adapter trimmed count map
        size_t* mCorrectionMatrix;                                 ///< 1x64 array to store various base to base correction accumulation numbers
        bool mSummarized;                                          ///< whether a FilterResult has been mSummarized(i.e mCorrectedBases calculated)
    public:
        /** Construct FilterResult object to store/calculate filter results
         * @param opt Options object including various filter options
         * @param paired if true the FilterResult object store pairend filter results
         */ 
        FilterResult(Options* opt, bool paired = false);
        
        /** Destruct FilterResult object, free mCorrectionMatrix */
        ~FilterResult();

        /** get mFilterReadStats array(reads filter status accumulation array)
         * @return pointer to mFilterReadStats
         */
        inline size_t* getFilterReadStats(){
            return mFilterReadStats;
        }

        /** add filter status result to filter status result accumulation array
         * @param result filter status
         */
        void addFilterResult(int result);

        /** add filter status result to filter status result accumulation array
         * @param result filter status
         * @param n reads number
         */
        void addFilterResult(int result, int n);

        /** add number of merged read pair number to filter results
         * @param n number of pairs merged
         */
        void addMergedPairs(int n);

        /** merge list of FilterResult object into one(combine filter results of various parts of a fastq file)
         * @param list a vector of pointer to FilterResult object
         * @return pointer to the merged FilterResult object
         */
        static FilterResult* merge(std::vector<FilterResult*>& list);
        
        /** Operator to output FilterResult to ostream
         * @param os reference of std::ostream object
         * @param re pointer of FilterResult object
         */
        friend std::ostream& operator<<(std::ostream& os, FilterResult* re);
        
        /** summary the FilterResult object, in case recalculate some parameters again and again
         * @param force if true, recalculate parameters
         */
        void summary(bool force = false);

        /** Update FilterResult when a adapter trimming happened
         * @param adapter adapter sequence trimmed
         * @param isR2 if true the adapter was trimmed from read2
         */
        void addAdapterTrimmed(const std::string& adapter, bool isR2 = false);

        /** Update FilterResult when adapter trimming happpened in pair end reads filtering
         * @param adapter1 adapter sequence trimmed from read1
         * @param adapter2 adapter sequence trimmed from read2
         */
        void addAdapterTrimmed(const std::string& adapter1, const std::string& adapter2);

        /** Report basic filter results to json file
         * @param j reference of json object
         */
        void reportJsonBasic(jsn::json& j);

        /** Report basic filter results to html file
         * @param ofs output file stream
         * @param totalReads total reads of fastq
         * @param totalBases total bases of fastq
         */
        void reportHtmlBasic(std::ofstream& ofs, size_t totalReads, size_t totalBases);

        /** Report detailed trimmed adapter sequence and counts in json fomat
         * @param j reference of json object
         * @param adapterCounts adapter sequence counts map
         */
        void reportAdaptersJsonDetails(jsn::json& j, std::map<std::string, size_t>& adapterCounts);
        
        /** Report summary trimmed adapter sequence and counts in json format
         * @param j reference of json object
         */
        void reportAdaptersJsonSummary(jsn::json& j);
        
        /** Report detailed trimmed adapter sequence and counts in html format
         * @param ofs output file stream
         * @param adapterCounts adapter sequence counts map
         * @param totalBases totalBases total bases of fastq
         */
        void reportAdaptersHtmlDetails(std::ofstream& ofs, std::map<std::string, size_t>& adapterCounts, size_t totalBases);

        /** Report summary trimmed adapter sequence and counts in html format
         * @param ofs output file stream
         * @param totalBases totalBases total bases of fastq
         */
        void reportAdaptersHtmlSummary(std::ofstream& ofs, size_t totalBases);

        /** get correction matrix(1x64 array)
         * @return pointer to mCorrectionMatrix
         */
        size_t* getCorrectionMatrix(){
            return mCorrectionMatrix;
        }

        /** get number of bases corrected
         * @return number of based corrected
         */
        size_t getTotalCorrectedBases();

        /** add one case of base correction to FilterResult object, \\nupdate the count in mCorrectionMatrix,\\n
         * mCorrectionMatrix[(from & 0x07) * 8 + (to & 0x07)] += 1
         * @param from source uncorrected base
         * @param to destination corrected base
         */
        void addCorrection(char from, char to);
       
        /** get base correction numbers in the pattern (from->to)
         * @param from source uncorrected base
         * @param to destination corrected base
         * @return base correction numbers in the pattern(from->to), just mCorrectionMatrix[(from & 0x07) * 8 + (to & 0x07)]
         */ 
        size_t getCorrectionNum(char from, char to);

        /** increase number of corrected reads
         * @param count number of corrected reads to increase
         */ 
        void incCorrectedReads(int count);

        /** add polyx trimmed statistics into results
         * @param base which base trimmed 0->A, 1->T, 2->C, 3->G
         * @param length number of base trimmmed
         */
        void addPolyXTrimmed(int base, int length);

        /** report polyx trimming statistics in json
         * @param j reference of json object
         */
        void reportPolyXTrimJson(jsn::json& j);
};

#endif
