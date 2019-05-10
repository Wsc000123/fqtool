#ifndef STATS_H
#define STATS_H

#include <map>
#include <iomanip>
#include <memory>
#include <vector>
#include <cstdio>
#include <string>
#include <cstdlib>
#include <sstream>
#include <ctml.hpp>
#include "read.h"
#include "util.h"
#include "options.h"
#include "evaluator.h"

/** Class to do statistics of a fastq file */
class Stats{
    private:
        Options* mOptions;                              ///< pointer to Options (including kmer and ora option objects)
        size_t mReads;                                  ///< total reads
        size_t mBases;                                  ///< total bases
        bool mIsRead2;                                  ///< this read in statistics is read2 if true
        int mMinReadLen;                                ///< minimal read length
        int mMaxReadLen;                                ///< maximum read length
        int mMinQual;                                   ///< minimum base quality
        int mMaxQual;                                   ///< maximum base quality
        int mEvaluatedSeqLen;                           ///< estimated read length
        int mCycles;                                    ///< maximum cycle
        int mBufLen;                                    ///< buffer length to store a static item of one read
        size_t mQ20Total;                               ///< number of bases with quality greater than 20
        size_t mQ30Total;                               ///< number of bases with quality greater than 30
        bool mSummarized;                               ///< summarized or not, if summarize() called, then mSummarized = true
        size_t mKmerMax;                                ///< the maximum mKmerLen kmer
        size_t mKmerMin;                                ///< the minimum mKmerLen kmer
        int mKmerLen = 5;                               ///< kmer length to calculate, default 5
        int mKmerBufLen;                                ///< kmer statistic array length, just equal 2 << (2 * mKmerLen)
        size_t mLengthSum;                              ///< total length of reads
        size_t mQ20Bases[8];                            ///< each base(ATCG) counts with quality greater than 20
        size_t mQ30Bases[8];                            ///< each base(ATCG) counts with quality greater than 30
        size_t mBaseContents[8];                        ///< each base(ATCG) counts
        size_t *mCycleQ20Bases[8];                      ///< each base(ATCG) counts at each cycle with quality greater than 20
        size_t *mCycleQ30Bases[8];                      ///< each base(ATCG) counts at each cycle with quality greater than 30
        size_t *mCycleBaseContents[8];                  ///< each base(ATCG) counts at each cycle
        size_t *mCycleBaseQuality[8];                   ///< each base(ATCG) quality at each cycle
        size_t *mCycleTotalBase;                        ///< total base counts of each cycle
        size_t *mCycleTotalQuality;                     ///< total quality of each cycle
        size_t *mKmer;                                  ///< kmer count array
        int mOverRepSampling = 100;                     ///< over representation analysis sampling frequence, default 100
        std::map<std::string, size_t> mOverReqSeqCount; ///< map of <repSeq, repSeqCount>, shoulde be initialized by Evaluator before statRead
        std::map<std::string, double*> mQualityCurves; ///< map of <statName, statNumber*> statNumber is a pointer to array of statistics of each cycle quality 
        std::map<std::string, double*> mContentCurves; ///< map of <statName, statNumber*> statNumber is a pointer to array of statistics of each cycle content
        std::map<std::string, size_t*> mOverRepSeqDist;///< map of <repseq, dist*>
    
    public:
        /** Construct a Stats object of fq 
         * @param opt pointer to Options object
         * @param isRead2 this fq is read2 if true
         * @param bufferMargin additional buffer length 
         */
        Stats(Options* opt, bool isRead2 = false, int bufferMargin = 1024);
        
        /** Destroy a Stats object Free allocated mamory */
        ~Stats();
        
        /** Allocate resources for a Stats object */
        void allocateRes();
        
        /** Initialize overrepresentation analysis
         */
        void initOverRepSeq();

        /** get minimum read length
         * @return minimum read length
         */
        int getMinReadLength();

        /** get maximum read length
         * @return maximum read length
         */
        int getMaxReadLength();

        /** get minimum base quality
         * @return minimum base quality
         */
        int getMinBaseQual();

        /** get maximum base quality
         * @return maximum base quality
         */
        int getMaxBaseQual();
         
        /** get cycle number
         * @return cycle number 
         */
        int getCycles();
        
        /** get total raeds number 
         * @return total read number 
         */
        size_t getReads();
        
        /** get total base number 
         * @return total base number 
         */
        size_t getBases();
        
        /** get number of bases with quality greater than 20 
         * @return number of bases with quality greater than 20
         */
        size_t getQ20();

        /** get number of bases with quality greater than 30  
         * @return number of bases with quality greater than 30 
         */
        size_t getQ30();

        /** get number of base G and C 
         * @return number of base G and C 
         */
        size_t getGCNumber();

        /** do statistics of one read 
         * @param r pointer to object Read
         */
        void statRead(Read* r);

        /** merge a list of Stats objects into one
         * @param list a list of Stats objects
         * @return a merged Stats object
         */
        static Stats* merge(std::vector<Stats*>& list);
        
        /** output Stats object to std::ostream
         * @param os std::ostream object
         * @param s reference of Stats object
         * @return std::ostream object
         */
        friend std::ostream& operator<<(std::ostream& os, const Stats& s);

        /** output Stats object to std::ostream
         * @param os std::ostream object
         * @param s pointer to Stats object
         * @return std::ostream object
         */
        friend std::ostream& operator<<(std::ostream& os, const Stats* s){
            os << (*s);
            return os;
        }
        
        /** Summary statistic items, get mBases, mCycles, mQ20Bases, mQ30Bases, mQ20Total, mQ30Total,
         * mQualityCurves["mean"], mQualityCurves["A/T/C/G/N"], mContentCurves["A/T/C/G/N"], mContentCurves["GC"]
         * mKmerMin, mKmerMax, after summary, summarize will be set to true
         * @param forced forced to run summarize
         * summarize will only run if (mSummarized is false) || (mSummarized is true && force = false)
         */
        void summarize(bool forced = false);
        
        /** Generate Html report
         * @param filteringType filtering type of report
         * @param readName library name 
         * @return a vector of nodes
         */
        std::vector<CTML::Node> reportHtml(std::string filteringType, std::string readName);

        /** Generate Quality part of the Html report
         * @param filteringType filtering type of report
         * @param readName library name 
         * @return quality node
         */
        CTML::Node reportHtmlQuality(std::string filteringType, std::string readName);
        
        /** Generate Content part of the Html report
         * @param filteringType filtering type of report
         * @param readName library name
         * @return contents node 
         */
        CTML::Node reportHtmlContents(std::string filteringType, std::string readName);
        
        /** Generate Kmer part of the Html report
         * @param filteringType filtering type of report
         * @param readName library name 
         * @return kmer node
         */
        CTML::Node reportHtmlKmer(std::string filteringType, std::string readName);
         
        /** Generate Over representation analysis  part of the Html report
         * @param filteringType filtering type of report
         * @param readName library name 
         * @return ora node
         */
        CTML::Node reportHtmlORA(std::string filteringType, std::string readName);
        
        /** whether is long read fq
         * @return true if cycles > 300
         */
        bool isLongRead();
        
        /** get mean read length 
         * @return mean read length 
         */
        int getMeanLength();

        /** convert an array of value type T to a string seperated by ","
         * @param list pointer to a T value array
         * @param size the length of the T value array
         */ 
        template<typename T>
        static std::string list2string(T* list, int size){
            std::stringstream ss;
            for(int i = 0; i < size; ++i){
                ss << list[i];
                if(i < size - 1){
                    ss << ",";
                }
            }
            return ss.str();
        }
       
        /** convert an array of T values to a string seperated by ","
         * @param list pointer to a T value array
         * @param size the length of the T value array
         * @param coords coordinates to define the ith value output
         * ith value output equals average of(list[coords[i-1]] ... list[coords[i]])
         */ 
        template<typename T>
        static std::string list2string(T* list, int size, size_t* coords){
            std::stringstream ss;
            long start = 0;
            long end = 0;
            T total = 0;
            for(int i = 0; i < size; ++i){
                if(i > 0){
                    start = coords[i - 1];
                }
                end = coords[i];
                total = 0;
                for(int j = start; j < end; ++j){
                    total += list[j];
                }
                if(end == start){
                    ss << "0";
                }else{
                    ss << total /(end - start);
                }
                if(i < size - 1){
                    ss << ",";
                }
            }
            return ss.str();
        }

    private:
        /** extend the array buffer for statistics longer
         * @param newBufLen the expected smallest buffer length
         * for performance, buffer length will increase to std::max(newBufLen + 100, 1.5 * newBufLen)
         */
        void extendBuffer(int newBufLen);
       
        /** make html popup value of kmer table
         * @param n the integer representation of kmer
         * @return popup value of kmer table
         */
        std::string makeKmerTD(size_t n);
        
        /** delete OverRepDist recources */ 
        void deleteOverRepSeqDist();

        /** test wheather count of seq is over represented
         * @param seq sequence to be evaluated
         * @param count seq count in sampling 
         * @return true if the seq is over represented
         */
        bool overRepPassed(const std::string& seq, size_t count);
};

#endif
