#ifndef WRITER_THREAD_H
#define WRITER_THREAD_H

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <mutex>
#include <atomic>
#include "util.h"
#include "writer.h"
#include "options.h"

/** class to hold a writer thread to write to one file from ring buffer */
class WriterThread{
    public:
        /** construct a WriterThread object
         * @param opt pointer to Options object
         * @param filename output filename
         */
        WriterThread(Options* opt, const std::string& filename);
       
        /** destroy a WriterThread object
         */ 
        ~WriterThread();
        
        /** construct a Writer object with output filename 
        * @param filename output filename of thie WriterThread
        */ 
        void iniWriter(const std::string& filename);
        
        /** construct a Writer object with output ofstream
         * @param ofs pointer to ofstream
         */
        void iniWriter(std::ofstream* ofs);
        
        /** construct a Writer object with a gzip file handler
         * @param gzfile a gZfile handler
         */
        void iniWriter(gzFile gzfile);

        /** free resources used by ringbuffer
         */ 
        void cleanup();

        /** test wheather this writing thread compoleted its work
         * @return true if mInputCompleted && (mOutputCounter == mInputCounter)
         */
        bool isCompleted();

        /** write mRingBuffer[mOutputCounter] to output if possible\n
         * if(if(mOutputCounter >= mInputCounter){waiting for 100s for input}\n
         * then try to write mRingBuffer[mOutputCounter] to output while(mOutputCounter < mInputCounter)\n
         * increase mOutputCounter after each writting
         */
        void output();

        /** feed C string to mRingBuffer and increase mInputCounter
         * @param cstr C string to add into mRingBuffer
         * @param size C string length
         */
        void input(char* cstr, size_t size);
        
        /** set mInputCompleted = true;
         * return true
         */
        bool setInputCompleted();
        
        /** get number of C strings in this thread ringbuffer to be written to output
         * @return mInputCounter - mOutputCounter
         */ 
        size_t bufferLength();
        
        /** get output fiilename
         * @return filename
         */
        inline std::string getFilename(){
            return mFilename;
        }

    private:
        /** destroy Writer object (Writer will close file if needed)
         */
        void deleteWriter();

    private:
        Options* mOptions;                  ///< pointer to Options
        Writer* mWriter;                    ///< Writer object to write cstring in ringbuffer into output
        std::string mFilename;              ///< output filename of this thread
        // for spliting control
        bool mInputCompleted;               ///< if true this thread have finished all writting
        std::atomic<size_t> mInputCounter;  ///< atomic type long integer to mark the index of input position in ringbuffer
        std::atomic<size_t> mOutputCounter; ///< atomic type long integer to mark the index of output position in ringbuffer
        char** mRingBuffer;                 ///< array to store pointer to C string to be written
        size_t* mRingBufferSizes;           ///< array to store the length of each C string to be written
};

#endif
