#include "writerthread.h"

WriterThread::WriterThread(Options* opt, const std::string& filename){
    mOptions = opt;
    mWriter = NULL;
    mInputCounter = 0;
    mOutputCounter = 0;
    mInputCompleted = false;
    mFilename = filename;

    mRingBuffer = new char*[mOptions->bufSize.maxPacksInReadPackRepo];
    std::memset(mRingBuffer, 0, sizeof(char*) * mOptions->bufSize.maxPacksInReadPackRepo);
    mRingBufferSizes = new size_t[mOptions->bufSize.maxPacksInReadPackRepo];
    std::memset(mRingBufferSizes, 0, sizeof(size_t) * mOptions->bufSize.maxPacksInReadPackRepo);
    iniWriter(mFilename);
}

WriterThread::~WriterThread(){
    cleanup();
    delete mRingBuffer;
    delete mRingBufferSizes;
}

bool WriterThread::isCompleted(){
    return mInputCompleted && (mOutputCounter == mInputCounter);
}

bool WriterThread::setInputCompleted(){
    mInputCompleted = true;
    return true;
}

void WriterThread::output(){
    util::loginfo("output started, start wait, mOutputCounter = " + std::to_string(mOutputCounter) + " mInputCounter = " + std::to_string(mInputCounter), mOptions->logmtx);
    if(mOutputCounter >= mInputCounter){
        usleep(100);
    }
    util::loginfo("output started, passed wait, mOutputCounter = " + std::to_string(mOutputCounter) + " mInputCounter = " + std::to_string(mInputCounter), mOptions->logmtx);
    while(mOutputCounter < mInputCounter){
        mWriter->write(mRingBuffer[mOutputCounter], 
                            mRingBufferSizes[mOutputCounter]);
        util::loginfo("write bytes: mRingBufferSizes[mOutputCounter] = " + std::to_string(mRingBufferSizes[mOutputCounter]), mOptions->logmtx);
        delete mRingBuffer[mOutputCounter];
        mRingBuffer[mOutputCounter] = NULL;
        ++mOutputCounter;
    }
    util::loginfo("output finished, mOutputCounter = " + std::to_string(mOutputCounter) + " mInputCounter = " + std::to_string(mInputCounter), mOptions->logmtx);
    mOutputCounter = mOutputCounter % mOptions->bufSize.maxPacksInReadPackRepo;
    mInputCounter = mInputCounter % mOptions->bufSize.maxPacksInReadPackRepo;
    util::loginfo("Counter Reseted, mOutputCounter = " + std::to_string(mOutputCounter) + " mInputCounter = " + std::to_string(mInputCounter), mOptions->logmtx);
}

void WriterThread::input(char* cstr, size_t size){
    while(mInputCounter == mOptions->bufSize.maxPacksInReadPackRepo){
        usleep(100);
    }
    if(mInputCounter < mOptions->bufSize.maxPacksInReadPackRepo){ 
        mRingBuffer[mInputCounter] = cstr;
        mRingBufferSizes[mInputCounter] = size;
        ++mInputCounter;
    }
}

void WriterThread::cleanup(){
    deleteWriter();
}

void WriterThread::deleteWriter(){
    if(mWriter != NULL){
        delete mWriter;
        mWriter = NULL;
    }
}

void WriterThread::iniWriter(const std::string& mFilename){
    deleteWriter();
    mWriter = new Writer(mFilename, mOptions->compression);
}

void WriterThread::iniWriter(std::ofstream* ofs){
    deleteWriter();
    mWriter = new Writer(ofs);
}

void WriterThread::iniWriter(gzFile gzfile){
    deleteWriter();
    mWriter = new Writer(gzfile);
}

size_t WriterThread::bufferLength(){
    return mInputCounter - mOutputCounter;
}
