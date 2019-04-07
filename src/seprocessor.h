#ifndef SE_PROCESSOR_H
#define SE_PROCESSOR_H

#include <string>
#include <memory>
#include <mutex>
#include <thread>
#include <atomic>
#include <unistd.h>
#include "util.h"
#include "read.h"
#include "polyx.h"
#include "stats.h"
#include "filter.h"
#include "fqreader.h"
#include "duplicate.h"
#include "umiprocessor.h"
#include "filterresult.h"
#include "writerthread.h"
#include "threadconfig.h"
#include "adaptertrimmer.h"


/** Struct to hold a bunch of redas pointers */
struct ReadPack{
    Read** data;   ///< array to store Read pointers
    int count;     ///< number of Read pointers in data
};

/** Struct to hold a bunch of ReadPack pointers */
struct ReadPackRepository{
    ReadPack** packBuffer;        ///< array to store ReadPack pointers
    std::atomic<size_t> readPos;    ///< the index in packBuffer next ReadPack pointer should read from
    std::atomic<size_t> writePos;   ///< the index in packBuffer next ReadPack pointer should put into
};

/** class to deal with single end fastq processing */
class SingleEndProcessor{
    public:
        /** construct a SingleEndProcessor object and set its default values 
         * @param opt pointer to Options
         */
        SingleEndProcessor(Options* opt);   
        
        /** destroy a SingleEndProcessor object and free resources */
        ~SingleEndProcessor();
       
        /** main control of single end fasq processing steps
         * @return true if a single end fastq process finished successfuly
         */ 
        bool process();
        
        /** process a pack of single end reads in one thread
         * @param pack pointer to a ReadPack object 
         * @param config a pointer to a ThreadConfig object 
         */
        void processSingleEnd(ReadPack* pack, ThreadConfig* config);
        
        /** initialize a ReadPackRepository object\n
         * make room for mRepo.packBuffer to store at most mOptions->bufSize.maxPacksInReadPackRepo packs in memory\n
         * and initialize mRepo.writePos = mRepo.readPos = 0;
         */
        void initReadPackRepository();

        /** destroy a ReadPackRepository and free allocated memory */
        void destroyReadPackRepository();
       
        /** add a fresh generated ReadPack pointer to ReadPackRepository
         * @param pack pointer to a ReadPack object
         */
        void producePack(ReadPack* pack);
       
        /** instruct a thread to consume a ReadPack\n
         * the thread will find the next proper mRepo.readPos\n
         * and extract ReacPack at mRepo.packBuffer[mRepo.readPos] to process\n
         * @param config pointer to ThreadConfig
         */
        void consumePack(ThreadConfig* config);
        
        /** a task(running asynchronously) to read fastq\n
         * fill a ReadPack and store the ReadPack into ReadPackRepository\n
         * continously till eof reached, but will pause at sometimes:\n
         * 1)mRepo.writePos - mRepo.readPos > mOptions->bufSize.maxPacksInMemory)
         * 2)readNum % (mOptions->bufSize.maxReadsInPack * mOptions->bufSize.maxPacksInMemory) == 0 && mLeftWriter)\n
         * && mLeftWriter->bufferLength() > mOptions->bufSize.maxPacksInMemory
         */ 
        void producerTask();
        
        /** a task(running asynchronously by each writing thread) to process some packs of Reads
         * @param config pointer to ThreadConfig pointer
         */
        void consumerTask(ThreadConfig* config);
        
        /** initialize ThreadConfig object for consumerTask running
         * @param config pointer to ThreadConfig object
         */
        void initConfig(ThreadConfig* config);
        
        /** will run only when split output is disabled\n
         * just create a WriterThread object for oubput writing and store it in mLeftWriter
         */
        void initOutput();
        
        /** will run only when split output is disabled\n
         * close WriterThread object used for output writing, just delete mLeftWriter
         */
        void closeOutput();
        
        /** will run only when split output is disabled\n
         * continously execute config->output() until finished\n
         * @param config pointer to a WriterThread
         */
        void writeTask(WriterThread* config);
        
        Options* mOptions;                   ///< pointer to Options
        ReadPackRepository mRepo;            ///< ReadPackRepository to store ReadPacks
        std::atomic<bool> mProduceFinished;  ///< if true the thread produce packs have finished work
        std::atomic<int> mFinishedThreads;   ///< number of threads who have finished their work
        std::mutex mInputMtx;                ///< mutex used to lock ReadPack extraction from ReadPackRepository by ThreadConfig
        std::mutex mOutputMtx;               ///< mutex used to lock WriterThread input when put one pack results into WriterThread 
        Filter* mFilter;                     ///< pointer to Filter to do various filter to each reads processed  
        gzFile mZipFile;                     ///< gzFile pointer used as output
        std::ofstream* mOutStream;           ///< filestream pointer used as output 
        UmiProcessor* mUmiProcessor;         ///< pointer to UmiProcessor to do umi processing
        WriterThread* mLeftWriter;           ///< pointer to WriterThread to perform writing if split output is disabled
        Duplicate* mDuplicate;               ///< pointer to Duplicate to do duplicate analysis
};

#endif
