#ifndef THREAD_CONFIG_H
#define THREAD_CONFIG_H

#include <string>
#include <vector>
#include "util.h"
#include "stats.h"
#include "writer.h"
#include "options.h"
#include "filterresult.h"

/** Class to configure a thread to write and store some results of the thread doing other than writing */
class ThreadConfig{
public:
    /** construct a ThreadConfig object
     * @param opt pointer to Options object
     * @param threadId this is just a manual created artifical thread marker, not the system allocated id
     * @param paired if true the input is paired end sequence fastq
     */
    ThreadConfig(Options* opt, int threadId, bool paired = false);
    
    /** destroy a ThreadConfig object
     */
    ~ThreadConfig();

    /** get the pointer to Stats object which are all statistic results of read 1 feeding in this thread before filtering
     * @return a pointer to Stats object
     */
    inline Stats* getPreStats1() {return mPreStats1;}

    /** get the pointer to Stats object which are all statistic results of read 1 feeding in this thread after filtering
     * @return a pointer to Stats object
     */
    inline Stats* getPostStats1() {return mPostStats1;}

    /** get the pointer to Stats object which are all statistic results of read 2 feeding in this thread before filtering
     * @return a pointer to Stats object
     */
    inline Stats* getPreStats2() {return mPreStats2;}
    
    /** get the pointer to Stats object which are all statistic results of read 2 feeding in this thread after filtering
     * @return a pointer to Stats object
     */
    inline Stats* getPostStats2() {return mPostStats2;}

    /** get the read1 Writer 
     * @return pointer to a Writer
     */
    inline Writer* getWriter1() {return mWriter1;}
    
    /** get the read2 Writer
     * @return pointer to a Writer
     */
    inline Writer* getWriter2() {return mWriter2;}
    
    /** get the FilterResult of this reads in this thread 
     * @return a pointer FilterResult object
     */
    inline FilterResult* getFilterResult() {return mFilterResult;}

    /** initialize a Writer with one filename
     * @param filename1 filename of read1
     */
    void initWriter(std::string filename1);
    
    /** initialize a Writer with two filenames
     * @param filename1 filename of read1
     * @param filename2 filename of read2
     */
    void initWriter(std::string filename1, std::string filename2);
    
    /** initialize a Writer with one ofstream
     * @param stream pointer to a ofstream
     */
    void initWriter(std::ofstream* stream);
    
    /** initialize a Writer with two ofstream
     * @param stream1 pointer to ofstream1
     * @param stream2 pointer to ofstream2
     */
    void initWriter(std::ofstream* stream1, std::ofstream* stream2);
    
    /** initialize a Writer with one gzFile
     * @param gzfile a gzFile handler
     */
    void initWriter(gzFile gzfile);
    
    /** initialize a Writer with two gzFile
     * @param gzfile1 a gzFile handler
     * @param gzfile2 a gzFile handler
     */
    void initWriter(gzFile gzfile1, gzFile gzfile2);
    
    /** add filter result of one read to be written in this thread into the FilterResult in this thread
     * @param result filter result returned by Filter
     */
    void addFilterResult(int result);

    /** get the manual crafted thread marker
     * @return mThreadId
     */
    int getThreadId() {return mThreadId;}
    
    // for splitting output
    /** increase the mCurrentSplitReads by readNum(written by this thread) and check if this thread should stop or continue
     * @param readNum number of reads this thread have written just now
     */
    void markProcessed(long readNum);
    
    /** initialize Writer for split output writting, just make new filenames containg the file part number\n
     * consist mOptions->split.digits 0 + "." + basename(mOptions->out1)
     */
    void initWriterForSplit();
    
    /** test whether this writing thread can be stoped now
     * mCurrentSplitReads >= mOptions->split.size && \n
     * !(mOptions->split.byFileLines || mWorkingSplit + mOptions->thread < mOptions->split.number) && \n
     * (mOptions->split.number % mOptions->thread >0 && mThreadId >= mOptions->split.number % mOptions->thread)
     * @return true if mCanBeStopped is true
     */
    bool canBeStopped();

    /** word do before destroy this object, just write empty files for spliting */
    void cleanup();

private:
    /** remove existing Writers and make them point to NULL */
    void deleteWriter();

    /** writing empty files for spliting */
    void writeEmptyFilesForSplitting();

private:
    Stats* mPreStats1;           ///< pointer to Stats object to hold prefilter read1 stats info
    Stats* mPostStats1;          ///< pointer to Stats object to hold afterfilter raed1 stats info
    Stats* mPreStats2;           ///< pointer to Stats object to hold prefilter read2 stats info
    Stats* mPostStats2;          ///< pointer to Stats object to hold afterfilter raed2 stats info
    Writer* mWriter1;            ///< pointer to read1 Writer object
    Writer* mWriter2;            ///< pointer to read2 Writer object
    Options* mOptions;           ///< pointer to Options object
    FilterResult* mFilterResult; ///< pointer to FilterResult
    // for spliting output
    int mThreadId;               ///< manual made artificial thread/split marker
    int mWorkingSplit;           ///< initial is just mThreadId, if this thread is full, it may increase mOptions->thread
    size_t mCurrentSplitReads;   ///< reads reads writen to current output part file
    bool mCanBeStopped;          ///< this thread have done work and can be stopped now if true
};

#endif
