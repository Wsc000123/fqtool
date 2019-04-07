#ifndef OVERLAP_ANALYSIS_H
#define OVERLAP_ANALYSIS_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "read.h"

/** Class to store overlap analysis results.*/
class OverlapResult{
    public:
        bool overlapped;  ///< true if read1/2 overlapped.
        int offset;       ///< offset value, 
                          ///< non-negative offset indicates overlap starts from 3'end of read2(TEMP_LEN > SEQ_LEN),\n
                          ///< negative offset indicates overlap starts from 5'end of read1(TEMP_LEN < SEQ_LEN).\n
        int overlapLen;   ///< overlap length of read1/2.
        int diff;         ///< number of different bases in the overlapped region of read1/2.
};

/** Class to do overlap analysis of Read1/2 or Seq1/2.*/
class OverlapAnalysis{
    public:
        /** Constructor of OverlapAnalysis */
        OverlapAnalysis() = default;

        /** Destructor of OverlapAnalysis */
        ~OverlapAnalysis() = default;

        /** Do overlap analysis of two Seq objects, seq1 and ~seq2 is overlapped if either one of the following two conditions satisfied\n
         * 1, overlap region >= overlapRequire && mismatch in overlap region <= overlapDiffLimit\n
         * 2, overlap region > 50 && mismatch in overlap region > overlapDiffLimit\n
         * @param s1 Seq object 
         * @param s2 Seq object
         * @param overlapDiffLimit maximum base differences allowed in the overlapped region
         * @param overlapRequire minimum required length of the overlapped region
         */
        static OverlapResult analyze(Seq& s1, Seq& s2, int overlapDiffLimit = 5, int overlapRequire = 30);

        /** Do overlap analysis of two Read objects, r1 and r2 is overlapped if either one of the following two conditions satisfied\n
         * 1, overlap region >= overlapRequire && mismatch in overlap region <= overlapDiffLimit\n
         * 2, overlap region > 50 && mismatch in overlap region > overlapDiffLimit\n
         * @param r1 pointer to Read object 
         * @param r2 pointer to Read object
         * @param overlapDiffLimit maximum base differences allowed in the overlapped region
         * @param overlapRequire minimum required length of the overlapped region
         */
        static OverlapResult analyze(Read* r1, Read* r2, int overlapDiffLimit = 5, int overlapRequire = 30);
};

#endif
