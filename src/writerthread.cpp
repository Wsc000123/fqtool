#include "writerthread.h"

WriterThread::WriterThread(Options* opt, const std::string& filename){
    mOptions = opt;
    mWriter = NULL;
    mInputCounter = 0;
    mOutputCounter = 0;
    mInputCompleted = false;
    mFilename = filename;

    mRingBuffer = new char*[COMMONCONST::MAX_PACKS_IN_READPACKREPO];
    std::memset(mRingBuffer, 0, sizeof(char*) * COMMONCONST::MAX_PACKS_IN_READPACKREPO);
    mRingBufferSizes = new size_t[COMMONCONST::MAX_PACKS_IN_READPACKREPO];
    std::memset(mRingBufferSizes, 0, sizeof(size_t) * COMMONCONST::MAX_PACKS_IN_READPACKREPO);
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
    if(mOutputCounter >= mInputCounter){
        usleep(100);
    }
    while(mOutputCounter < mInputCounter){
        mWriter->write(mRingBuffer[mOutputCounter], 
                            mRingBufferSizes[mOutputCounter]);
        delete mRingBuffer[mOutputCounter];
        mRingBuffer[mOutputCounter] = NULL;
        ++mOutputCounter;
    }
    mOutputCounter = mOutputCounter % COMMONCONST::MAX_PACKS_IN_READPACKREPO;
    mInputCounter = mInputCounter % COMMONCONST::MAX_PACKS_IN_READPACKREPO;
}

void WriterThread::input(char* cstr, size_t size){
    while(mInputCounter == COMMONCONST::MAX_PACKS_IN_READPACKREPO){
        usleep(100);
    }
    if(mInputCounter < COMMONCONST::MAX_PACKS_IN_READPACKREPO){ 
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
