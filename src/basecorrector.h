#ifndef BASE_CORRECTOR_H
#define BASE_CORRECTOR_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include "util.h"
#include "overlapanalysis.h"
#include "filterresult.h"

/** Class to do base correction */
class BaseCorrector{
    public:
        /** Construct a BaseCorrector object */
        BaseCorrector();
        /** Destroy a BaseCorrector object */
        ~BaseCorrector();

        /** correct pair of reads by overlap analysis\n
         * will use overlap analysis result to correct bases if only at most 5 mismatches happened\n
         * @param r1 pointer to Read object(read1)
         * @param r2 pointer to Read object(read2)
         * @param fr pointer to FilterResult object(may be NULL)
         */
        static int correctByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr);
        
        /** correct pair of reads by overlap analysis\n
         * will use overlap analysis result to correct bases if only at most 5 mismatches happened\n
         * @param r1 pointer to Read object(read1)
         * @param r2 pointer to Read object(read2)
         * @param fr pointer to FilterResult object(may be NULL)
         * @param ov reference of OverlapResult object
         */
        static int correctByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr, OverlapResult& ov);
};

#endif
