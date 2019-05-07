#ifndef DUPLICATE_H
#define DUPLICATE_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <memory>
#include <mutex>
#include <cmath>
#include "options.h"
#include "read.h"
#include "overlapanalysis.h"

/** Class to do reads duplication analysis */
class Duplicate{
    public:
        /** Construct a Duplicate object and allocate resources */
        Duplicate(Options* opt);
        /** Remove a Duplicate object and free resources */
        ~Duplicate();

        /** Calculate key, kmer32 and gc of one read\n
         * calculate key, first mKeyLenInBase sequence to uint64_t and truncated to uint32_t\n
         * calculate kmer32, (readLength - 37) to (readLength - 5 sequence) || first 32 bases to uint64_t\n
         * calculate gc, gc ratio of whole read * 255.0\n
         * update stats array
         */
        void statRead(Read* r1);

        /** Calculate key, kmer32 and gc of a pair of read\n
         * calculate key use first mKeyLenInBase sequence of read1\n
         * calculate kmer32 use first 32 bases of read2\n
         * calculate gc us both read1/read2\n
         * update stats array
         */
        void statPair(Read* r1, Read* r2);

        /** Add duplicate statistical items to record array\n
         * mDups[key] = (the kmer32 to key with the highest GC encountered)\n
         * mCounts[key] = (the number of kmer32 with the same key)\n
         * mGC[key] = (the GC ratio * 255.0 of the kmer32 to key)\n
         */
        void addRecord(uint32_t key, uint64_t kmer32, uint8_t gc);
        
        /** Do analysis of all reads\n
         * @param hist hist[i] is corresponding the mCounts value i, if greater than histSize just keep it in histSize - 1
         * @param meanGC meanGC[i] is the corresponding average gc ratio of all keys with hist[i] mCounts value
         * @param histSize the length of hist/meanGC array 
         * @return total duplicate ratio of all reads
         */
        double statAll(size_t* hist, double* meanGC, size_t histSize);
        
        /** Convert an sequence consists of ATCG in to integer
         * @param cstr a C string pointer
         * @param start starting position of cstr to convert
         * @param keylen length of cstr substring to convert
         * @param valid if true the convert is valid, invalidate conversion only will happen if any N in cstr[start]...cstr[start + keylen - 1]
         */ 
        uint64_t seq2int(const char* cstr, int start, int keylen, bool& valid);
    
    private: 
        Options* mOptions;     ///< Options Object to provide duplicate analysis options
        int mKeyLenInBase;     ///< the length of the key in bases
        uint64_t mKeyLenInBit; ///< the bits number needed to represent all kinds of keys
        uint64_t* mDups;       ///< mDups[key] = (the kmer32 with the highest GC encountered)
        uint32_t* mCounts;     ///< mCounts[key] = (the same kmer32 number eccountered)
        uint8_t* mGC;          ///< mGC[key] = (the GC ratio * 255.0 of the kmer32)
        std::mutex mAddLock;   ///< lock used to add record
};
        
#endif
