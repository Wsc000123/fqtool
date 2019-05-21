#ifndef POLY_X_H
#define POLY_X_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include "read.h"
#include "filterresult.h"

/** Class to do Poly X trimming */
class PolyX{
    public:
        /** construct a PolyX object */
        PolyX();
        /** Destroy a PolyX object */
        ~PolyX();

        /** trim polyG from 3' end\n
         * if at most 5 mismatch in length l, and l >= compareReq, trim it
         * if at most 1 mismatch for each 8 bases, and l >= compareReq, trim it
         * @param r1 pointer to Read object(read1)
         * @param r2 pointer to Read object(read2)
         * @param compareReq required length of sequence to be polyG
         * @param maxMismatch max mismatches allowed for a sequence against G
         * @param allowedOneMismatchForEach max mismatches allowed for each allowedOneMismatchForEach bases against G
         * @param fr pointer to FilterResult object
         */
        static void trimPolyG(Read* r1, Read* r2, int compareReq, int maxMismatch, int allowedOneMismatchForEach, FilterResult* fr = NULL);
        
        /** trim polyG from 3' end\n
         * if at most 5 mismatch in length l, and l >= compareReq, trim it
         * if at most 1 mismatch for each 8 bases, and l >= compareReq, trim it
         * @param r pointer to Read object
         * @param compareReq required length of sequence to be polyG
         * @param maxMismatch max mismatches allowed for a sequence against G
         * @param allowedOneMismatchForEach max mismatches allowed for each allowedOneMismatchForEach bases against G
         * @param fr pointer to FilterResult object
         */
        static void trimPolyG(Read* r, int compareReq, int maxMismatch, int allowedOneMismatchForEach, FilterResult* fr = NULL);
        
        /** trim polyX from 3' end\n
         * if at most 5 mismatch in length l, and l >= compareReq, trim it
         * if at most 1 mismatch for each 8 bases, and l >= compareReq, trim it
         * @param r1 pointer to Read object(read1)
         * @param r2 pointer to Read object(read2)
         * @param trimChr nucleotides to be trimmed
         * @param compareReq required length of sequence to be polyX
         * @param maxMismatch max mismatches allowed for a sequence against X
         * @param allowedOneMismatchForEach max mismatches allowed for each allowedOneMismatchForEach bases against X
         * @param fr pointer to FilterResult object
         */
        static void trimPolyX(Read* r1, Read* r2, std::string trimChr,  int compareReq, int maxMismatch, int allowedOneMismatchForEach, FilterResult* fr = NULL);

        /** trim polyX from 3' end\n
         * if at most 5 mismatch in length l, and l >= compareReq, trim it
         * if at most 1 mismatch for each 8 bases, and l >= compareReq, trim it
         * @param r pointer to Read object
         * @param trimChr nucleotides to be trimmed
         * @param compareReq required length of sequence to be polyX
         * @param maxMismatch max mismatches allowed for a sequence against X
         * @param allowedOneMismatchForEach max mismatches allowed for each allowedOneMismatchForEach bases against X
         * @param fr pointer to FilterResult object
         */
        static void trimPolyX(Read* r, std::string trimChr, int compareReq, int maxMismatch, int allowedOneMismatchForEach, FilterResult* fr = NULL);
};

#endif
