#ifndef FQREADER_H
#define FQREADER_H

#include <zlib.h>
#include "util.h"
#include "read.h"
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <iostream>

/** Class to hold a fastq reader
 */
class FqReader{
    std::string mFileName;  ///< name of fastq file
    gzFile mGzipFile;       ///< gzFile to store opened mZipped file handler
    FILE* mFile;            ///< FILE pointer to store opened plain file handler
    bool mZipped;           ///< the fastq file is gzipped if true
    bool mHasQuality;       ///< the fastq file has quality sequence if true
    bool mPhread64;         ///< the fastq quality sequence is encoded as ASCII 64 based if true
    char* mBuf;             ///< mBuffer to store a chunk of characters read from mGzipFile or mFile
    int mBufDataLen;        ///< the length of characters read into the mBuffer after the last read
    int mBufUsedLen;        ///< the length of characters already consumed in the mBuffer 
    bool mStdinMode;        ///< read from stdin if true
    bool mNoLineBreakAtEnd; ///< the fastq file has no '\n' as a line break at the last line if true
    int mFqBufSize;         ///< the mBuffer size used to read
    
    public:
        /** Construct a FqReader with filename and hasQuality, phread64 arguments
         * @param filename Name of the fastq file
         * @param hasQuality The fastq has quality sequence if true
         * @param phread64 The fastq quality sequence is encoded as ASCII 64 if true
         */
        FqReader(const std::string& filename, const bool& hasQuality = true, const bool& phread64 = false);
        
        /** FqReader Destructor
         */
        ~FqReader();
       
        /** Tell whether the fastq file is zipped or not
         * @return true if the fasq file is zipped 
         */
        bool isZipped();

        /** Get the number of Bytes read from(write to) the fastq file in the meantime, store it in bytesRead
         *  Get the total number of Bytes in the fastq file, store it in bytesTotal
         */
        void getBytes(size_t& bytesRead, size_t& bytesTotal);

        /** Try to get the next Read record from the current FqReader
         * @return a pointer to object Read read from the current FqReader or NULL if failed during reading
         */
        Read* read();

        /** Tell whether the FqReader has reach the endof file
         * @return true if eof reached
         */
        bool eof();

        /** Tell whether the fastq file has no '\n' as a line break at the last line
         * @return true if the fastq file has no '\n' as a line break at the last line
         */ 
        bool hasNoLineBreakAtEnd();

    public:

        /** tell whether the file is a mZipped fastq file
         * @param filename file name
         * @return true if the file is a mZipped fastq file(with suffix ".fastq.gz", ".fasta.gz", ".fq.gz" or ".fa.gz")
         */
        static bool isZippedFq(const std::string& filename);
        
        /** tell whether the file is a fastq format file
         * @param filename file name
         * @return true if the file is a fastq file(with suffix ".fastq", ".fasta", ".fq" or ".fa)
         */
        static bool isfq(const std::string& filename);
    
        /** get just one line from the mBuf, update the mBuf if needed
         */
        std::string getLine();

    private:

        /** initialize the FqReader:
         * 1, open file and store file handler into mGzipFile or mFile or read from stdin
         * 2, set the starting position for the next read on compressed file stream file to the beginning of file 
         * 3, update the file format mZipped
         * 4, call readToBuf() to try to fill the mBuf from first reading
         */ 
        void init();
        
        /** close the FqReader:
         * 1, close file handler
         * 2, set file handler to NULL
         */
        void close();

        /** trim \n, \r or \r\n in the tail of the line
         */
        void clearLineBreaks(char* line);
        
        /** try a reading action to fill the mBuf,
         * update mBufDataLen to Bytes read during this reading action
         * reset mBufUsedLen to zero
         * if read the last line(mBuf is not filled), update mNoLineBreakAtEnd
         */
        void readToBuf();
};

/** Class to hold a pair end fastq reader
 */
class FqReaderPair{
    public:
        FqReader* left;   ///< Pointer to FqReader of read1
        FqReader* right;  ///< Pointer to FqReader of read2
        bool mInterleaved; ///< Whether the pairend fastq is interleaved(r1/r2 combined into one file) or not(r1/r2 in seperated files)

    public:
        /** Construct a FqReaderPair from two pointers of FqReader
         * @param rl Pointer to FqReader of read1
         * @param rr Pointer to FqReader of read2
         */
        FqReaderPair(FqReader* rl, FqReader* rr) : left(rl), right(rr){}
        
        /** Construct a FqReaderPair from file names
         * @param lname read1 filename
         * @param rname read2 filename
         * @param hasQual whether the fastq has quality sequence or not, default true
         * @param mPhread64 whether the fastq quality sequence is encoded in ASCII 64 based chacacters, default false
         * @param interleaved whether the pairend fastq is interleaved(r1/r2 combined into one file) or not(r1/r2 in seperated files), default false
         */
        FqReaderPair(const std::string& lname, const std::string& rname, const bool& hasQual = true,
                     const bool& mPhread64 = false, const bool& interleaved = false);
        
        /** Destructor of FqReaderPair
         */ 
        ~FqReaderPair();

        /** try to read to get a pair of read from fastq
         * @return Pointer to object of ReadPair
         */
        ReadPair* read();
};

#endif
