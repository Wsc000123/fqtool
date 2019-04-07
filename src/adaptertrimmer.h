#ifndef ADAPTER_TRIMMER_H
#define ADAPTER_TRIMMER_H

#include <cstdlib>
#include <cstdio>
#include <string>
#include "overlapanalysis.h"
#include "filterresult.h"

/** Class to do adapter trimming */
class AdapterTrimmer{
    public:
        /** Construct an AdapterTrimmer */
        AdapterTrimmer();
        /** Destroy an AdapterTrimmer */
        ~AdapterTrimmer();

        /** Trim adapter by overlap analysis (auao detect adapter sequence)
         * @param r1 pointer to Read object(read1)
         * @param r2 pointer to Read object(read2)
         * @param fr pointer to FilterResult object
         */
        static bool trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr);
        
        /** Trim adapter by overlap analysis (based on overlap analysis results)\n
         * maximum mismatch is 5 and minimum overlap length is 1/3 the read length
         * @param r1 pointer to Read object(read1)
         * @param r2 pointer to Read object(read2)
         * @param fr pointer to FilterResult object
         * @param ov referenct to OverlapResult
         */ 
        static bool trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr, OverlapResult& ov);

        /** Trim adapter by provided adapter sequence\n
         * for every 8 bases overlap at most 1 mismatch is allowed
         * @param r pointer to Read object
         * @param fr pointer to FilterResult object
         * @param adapter adapter sequence provided externally
         * @param isR2 this is read2 of a pe fq if true
         */
        static bool trimBySequence(Read* r, FilterResult* fr, std::string& adapter, bool isR2 = false);
};

#endif
