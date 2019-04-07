#ifndef PE_PROCESSOR_H
#define PE_PROCESSOR_H

#include <string>
#include <memory>
#include <mutex>
#include <thread>
#include <atomic>
#include <unistd.h>
#include "read.h"
#include "util.h"
#include "filter.h"
#include "polyx.h"
#include "stats.h"
#include "common.h"
#include "fqreader.h"
#include "duplicate.h"
#include "umiprocessor.h"
#include "writerthread.h"
#include "threadconfig.h"
#include "filterresult.h"
#include "basecorrector.h"
#include "adaptertrimmer.h"

/** struct to store pointers of ReadPair */
struct ReadPairPack {
    ReadPair** data; ///< array to store ReadPair pointers
    int count;       ///< number of RaedPair pointers stored in data
};
   
/** struct to store pointers of ReadPairPack */
struct ReadPairPackRepository{
    ReadPairPack** packBuffer;  ///< array to store ReadPairPack pointers
    std::atomic<size_t> readPos;  ///< atomic type long integer to store the index next pointer of ReadPairPack will be extracted
    std::atomic<size_t> writePos; ///< atomic type long integer to store the index next pointer of ReadPairPack will be writen in
};

/** class to process pair end fastq */
class PairEndProcessor{
    public:
        /** construct a PairEndProcessor and set default values
         * @param opt pointer to Options
         */
        PairEndProcessor(Options* opt);
        
        /** destroy a PairEndProcessor object */
        ~PairEndProcessor();

        /** process a pair of fastq
         * @return true if process finished successfully
         */
        bool process();

        /** process a ReadPairPack 
         * @param pack pointer to ReadPairPack
         * @param config pointer to ThreadConfig
         * @return true if finish process
         */
        bool processPairEnd(ReadPairPack* pack, ThreadConfig* config);
        
        /** initialize ReadPairPackRepository\n
         * allocate memory to store at most mOptions->bufSize.maxPacksInReadPackRepo ReadPairPack pointers
         * set readPos and writePos of ReadPairPackRepository to be zero
         */
        void initReadPairPackRepository();
        
        /** destroy ReadPairPackRepository and free memory used */
        void destroyReadPairPackRepository();

        /** put a newly generated ReadPairPack pointer to the proper position in ReadPairPackRepository\n
         * that is writePos in packBuffer and increase writePos
         * @param pack pointer to ReadPairPack
         */
        void producePack(ReadPairPack* pack);

        /** extract a ReadPairPack from ReadPairPackRepository and process in a thread
         * @param config pointer to ThreadConfig
         */
        void consumePack(ThreadConfig* config);

        /** a task running asynchronously to read pair of reads into memory\n
         * put into a ReadPairPack firstly, if ReadPairPack filled or finished reading\n
         * put this ReadPairPack into the ReadPairPackRepository\n
         */
        void producerTask();

        /** a task running asynchronously to process read pairs in a thread
         * @param config pointer to ThreadConfig
         */
        void consumerTask(ThreadConfig* config);

        /** initialize a thread for process read pairs
         * @param config poiter to ThreadConcig
         */
        void initConfig(ThreadConfig* config);

        /** initialize two WriterThread for output if split output is disabled\n
         * and store pointer to this WriterThread into mLeftWriter and mRightWriter
         */
        void initOutput();

        /** close two WriterThread after writing if split output is disabled */
        void closeOutput();

        /** calculate insertsize of a pair of reads
         * @param r1 pointer to Read object (read1)
         * @param r2 pointer to Read object (read2)
         * @param ov reference to OverlapResult object
         */
        void statInsertSize(Read* r1, Read* r2, OverlapResult& ov);

        /** get insert size peak
         * @return inert size peak
         */
        int getPeakInsertSize();

        /** writing task running asynchronously to write pair of reads results to output
         * @param config pointer to a WriterThread object
         */
        void writeTask(WriterThread* config);

    private:
        Options* mOptions;                   ///< a pointer to object Options
        ReadPairPackRepository mRepo;        ///< ReadPairPackRepository object to store pointers of ReadPairPack
        std::atomic<bool> mProduceFinished;  ///< an atom type bool value to mark all reads have been read if true
        std::atomic<int> mFinishedThreads;   ///< an atom type int value to store the finished writing threads number
        std::mutex mOutputMtx;               ///< a mutex object to be locked when mRepo is extracted to be processed
        std::mutex mInputMtx;                ///< a mutex object to be locked when input to the WriterThread 
        Filter* mFilter;                     ///< a pointer to a Filter object to do various filter of pe reads
        gzFile mZipFile1;                    ///< gzFile to output read1 results
        gzFile mZipFile2;                    ///< gzFile to output read2 results
        std::ofstream* mOutStream1;          ///< output filestream to output read1
        std::ofstream* mOutStream2;          ///< output filestream to output read2
        UmiProcessor* mUmiProcessor;         ///< pointer to UmiProcessor to do umi process
        long* mInsertSizeHist;               ///< array to store different insert size counts
        WriterThread* mLeftWriter;           ///< pointer to a WriterThread object to write read1
        WriterThread* mRightWriter;          ///< pointer to a WriterThread object to write read2
        Duplicate* mDuplicate;               ///< pointer to a Duplicate object to du duplicate analysis
};

#endif
