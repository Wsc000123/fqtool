#include "seprocessor.h"

SingleEndProcessor::SingleEndProcessor(Options* opt){
    mOptions = opt;
    mProduceFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(mOptions);
    mOutStream = NULL;
    mZipFile = NULL;
    mUmiProcessor = new UmiProcessor(mOptions);
    mLeftWriter = NULL;
    mFailedWriter = NULL;
    mDuplicate = NULL;
    if(mOptions->duplicate.enabled){
        mDuplicate = new Duplicate(mOptions);
    }
}

SingleEndProcessor::~SingleEndProcessor(){
    delete mFilter;
    delete mUmiProcessor;
    if(mDuplicate){
        delete mDuplicate;
        mDuplicate = NULL;
    }
}

void SingleEndProcessor::initOutput(){
    if(!mOptions->failedOut.empty()){
        mFailedWriter = new WriterThread(mOptions, mOptions->failedOut);
    }
    if(mOptions->out1.empty()){
        return;
    }
    mLeftWriter = new WriterThread(mOptions, mOptions->out1);
}

void SingleEndProcessor::closeOutput(){
    if(mLeftWriter){
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
    if(mFailedWriter){
        delete mFailedWriter;
        mFailedWriter = NULL;
    }
}

void SingleEndProcessor::initConfig(ThreadConfig* config){
    if(mOptions->out1.empty()){
        return;
    }
    
    if(mOptions->split.enabled){
        config->initWriterForSplit();
    }
}

void SingleEndProcessor::initReadPackRepository(){
    mRepo.packBuffer = new ReadPack*[mOptions->bufSize.maxPacksInReadPackRepo];
    std::memset(mRepo.packBuffer, 0, sizeof(ReadPack*) * mOptions->bufSize.maxPacksInReadPackRepo);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
}

void SingleEndProcessor::destroyReadPackRepository(){
    delete mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void SingleEndProcessor::producePack(ReadPack* pack){
    while(mRepo.writePos ==  mOptions->bufSize.maxPacksInReadPackRepo){
        usleep(1);
    }
    mRepo.packBuffer[mRepo.writePos] = pack;
    util::loginfo("producer produced pack " + std::to_string(mRepo.writePos), mOptions->logmtx);
    ++mRepo.writePos;
}

void SingleEndProcessor::producerTask(){
    util::loginfo("loading data started", mOptions->logmtx);
    size_t readNum = 0;
    Read** data = new Read*[mOptions->bufSize.maxReadsInPack];
    std::memset(data, 0, sizeof(Read*) * mOptions->bufSize.maxReadsInPack);
    FqReader reader(mOptions->in1, true, mOptions->phred64);
    size_t count = 0;
    while(true){
        Read* read = reader.read();
        ++readNum;
        if(!read){
            if(count == 0){
                break;
            }
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            data = NULL;
            if(read){
                delete read;
                read = NULL;
            }
            mProduceFinished = true;
            break;
        }
        data[count] = read;
        ++count;
        if(count == mOptions->bufSize.maxReadsInPack){
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            data = new Read*[mOptions->bufSize.maxReadsInPack];
            std::memset(data, 0, sizeof(Read*) * mOptions->bufSize.maxReadsInPack);
            while(mRepo.writePos - mRepo.readPos > mOptions->bufSize.maxPacksInMemory){
                usleep(100);
            }
            readNum += count;
            if(readNum % (mOptions->bufSize.maxReadsInPack * mOptions->bufSize.maxPacksInMemory) == 0 && mLeftWriter){
                while(mLeftWriter->bufferLength() > mOptions->bufSize.maxPacksInMemory){
                    usleep(100);
                }
            }
            count = 0;
        }
    }
    mProduceFinished = true;
    util::loginfo("loaded reads: " + std::to_string(readNum - 1), mOptions->logmtx);
}

void SingleEndProcessor::consumerTask(ThreadConfig* config){
    while(true){
        if(config->canBeStopped()){
            ++mFinishedThreads;
            break;
        }
        while(mRepo.writePos <= mRepo.readPos){
            if(mProduceFinished){
                break;
            }
            usleep(1);
        }
        if(mProduceFinished && mRepo.writePos == mRepo.readPos){
            ++mFinishedThreads;
            break;
        }
        consumePack(config);
    }
    if(mFinishedThreads == mOptions->thread){
        if(mLeftWriter){
            mLeftWriter->setInputCompleted();
        }
        if(mFailedWriter){
            mFailedWriter->setInputCompleted();
        }
    }
    util::loginfo("thread " + std::to_string(config->getThreadId()) + " finished", mOptions->logmtx);
}

void SingleEndProcessor::consumePack(ThreadConfig* config){
    mInputMtx.lock();
    while(mRepo.writePos <= mRepo.readPos){
        usleep(1);
        if(mProduceFinished){
            mInputMtx.unlock();
            return;
        }
    }
    size_t packNum = mRepo.readPos;
    ReadPack* data = mRepo.packBuffer[packNum];
    ++mRepo.readPos;
    if(mRepo.readPos == mOptions->bufSize.maxPacksInReadPackRepo){
        mRepo.readPos = mRepo.readPos % mOptions->bufSize.maxPacksInReadPackRepo;
        mRepo.writePos = mRepo.writePos % mOptions->bufSize.maxPacksInReadPackRepo;
    }
    mInputMtx.unlock();
    util::loginfo("thread " + std::to_string(config->getThreadId()) + " start processing pack " + std::to_string(packNum), mOptions->logmtx);
    processSingleEnd(data, config);
    util::loginfo("thread " + std::to_string(config->getThreadId()) + " finish processing pack " + std::to_string(packNum), mOptions->logmtx);
}

bool SingleEndProcessor::process(){
    if(!mOptions->split.enabled){
        initOutput();
    }
    initReadPackRepository();
    util::loginfo("read pack repo initialized", mOptions->logmtx);
    std::thread producer(std::bind(&SingleEndProcessor::producerTask, this));
    util::loginfo("producer thread started", mOptions->logmtx);
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t = 0; t < mOptions->thread; ++t){
        configs[t] = new ThreadConfig(mOptions, t, false);
        initConfig(configs[t]);
    }
    std::thread** threads = new std::thread*[mOptions->thread];
    for(int t = 0; t < mOptions->thread; ++t){
        threads[t] = new std::thread(std::bind(&SingleEndProcessor::consumerTask, this, configs[t]));
    }
    util::loginfo(std::to_string(mOptions->thread) + " working threads started", mOptions->logmtx);
    std::thread* leftWriterThread = NULL;
    if(mLeftWriter){
        leftWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mLeftWriter));
        util::loginfo("read1 writer thread started", mOptions->logmtx);
    }
    std::thread* failedWriterThread = NULL;
    if(mFailedWriter){
        failedWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mFailedWriter));
        util::loginfo("failed reads writer thread started", mOptions->logmtx);
    }
    producer.join();
    util::loginfo("producer thread finished", mOptions->logmtx);
    for(int t = 0; t < mOptions->thread; ++t){
        threads[t]->join();
    }
    util::loginfo("working threads finished", mOptions->logmtx);
    if(!mOptions->split.enabled){
        if(leftWriterThread){
            leftWriterThread->join();
            util::loginfo("read1 writer thread finished", mOptions->logmtx);
        }
    }
    if(mFailedWriter){
        failedWriterThread->join();
        util::loginfo("failed reads writer thread finished", mOptions->logmtx);
    }
    util::loginfo("start generating reports", mOptions->logmtx);
    std::vector<Stats*> preStats;
    std::vector<Stats*> postStats;
    std::vector<FilterResult*> filterResults;
    for(int t = 0; t < mOptions->thread; ++t){
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
        filterResults.push_back(configs[t]->getFilterResult());
    }
    Stats* finalPreStats = Stats::merge(preStats);
    Stats* finalPostStats = Stats::merge(postStats);
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);
    // output filter results
    size_t* dupHist = NULL;
    double* dupMeanGC = NULL;
    double dupRate = 0.0;
    // output duplicate results
    if(mOptions->duplicate.enabled){
        dupHist = new size_t[mOptions->duplicate.histSize];
        std::memset(dupHist, 0, sizeof(int) * mOptions->duplicate.histSize);
        dupMeanGC = new double[mOptions->duplicate.histSize];
        std::memset(dupMeanGC, 0, sizeof(double) * mOptions->duplicate.histSize);
        dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
    }
    JsonReporter jr(mOptions);
    jr.setDupHist(dupHist, dupMeanGC, dupRate);
    jr.report(finalFilterResult, finalPreStats, finalPostStats);
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.report(finalFilterResult, finalPreStats, finalPostStats);
    util::loginfo("finish generating reports", mOptions->logmtx);
    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    delete finalPreStats;
    delete finalPostStats;
    delete finalFilterResult;

    if(mOptions->duplicate.enabled) {
        delete[] dupHist;
        delete[] dupMeanGC;
    }

    delete[] threads;
    delete[] configs;

    if(leftWriterThread){
        delete leftWriterThread;
    }
    if(failedWriterThread){
        delete failedWriterThread;
    }
    if(!mOptions->split.enabled){
        closeOutput();
    }

    return true;
}

void SingleEndProcessor::processSingleEnd(ReadPack* pack, ThreadConfig *config){
    std::string outstr;
    std::string failedOut;
    int readPassed = 0;
    for(int p = 0; p < pack->count; ++p){
        // original read1
        Read* or1 = pack->data[p];
        // stats the original read before trimming 
        config->getPreStats1()->statRead(or1);
        // handling the duplication profiling
        if(mOptions->duplicate.enabled){
            mDuplicate->statRead(or1);
        }
        // filter by index
        if(mOptions->indexFilter.enabled && mFilter->filterByIndex(or1)){
            delete or1;
            continue;
        }
        // umi processing
        if(mOptions->umi.enabled){
            mUmiProcessor->process(or1);
        }
        // trim in head and tail, and apply quality cut in sliding window
        Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1);
        // polyG trimming
        if(r1 != NULL){
            if(mOptions->polyGTrim.enabled){
                PolyX::trimPolyG(r1, mOptions->polyGTrim.minLen, mOptions->polyGTrim.maxMismatch, mOptions->polyGTrim.allowedOneMismatchForEach, config->getFilterResult());
            }
        }
        // adapter trimming
        if(r1 != NULL && mOptions->adapter.enableTriming && mOptions->adapter.adapterSeqR1Provided){
            AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.inputAdapterSeqR1);
        }
        // polyX trimming
        if(r1 != NULL){
            if(mOptions->polyXTrim.enabled){
                PolyX::trimPolyX(r1, mOptions->polyXTrim.trimChr, mOptions->polyXTrim.minLen, 
                                 mOptions->polyXTrim.maxMismatch, mOptions->polyXTrim.allowedOneMismatchForEach, config->getFilterResult());
            }
        }
        // trim max length
        if(r1 != NULL){
            if(mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length()){
                r1->resize(mOptions->trim.maxLen1);
            }
        }

        // get quality quality nbase length complexity...passing status
        int result = mFilter->passFilter(r1);
        config->addFilterResult(result);
        // stats the read after filtering
        if(r1 != NULL && result == COMMONCONST::PASS_FILTER){
            outstr += r1->toString();
            config->getPostStats1()->statRead(r1);
            ++readPassed;
        }else if(mFailedWriter){
            failedOut += or1->toStringWithTag(COMMONCONST::FAILED_TYPES[result]);
        }
        // cleanup memory
        delete or1;
        if(r1 != or1 && r1 != NULL){
            delete  r1;
        }
    }
    // if splitting output, then no lock is need since different threads write different files
    if(!mOptions->split.enabled){
        mOutputMtx.lock();
    }
    if(mOptions->outputToSTDOUT){
        std::fwrite(outstr.c_str(), 1, outstr.length(), stdout);
    }else if(mOptions->split.enabled){
        // split output by each worker thread
        if(!mOptions->out1.empty()){
            config->getWriter1()->writeString(outstr);
        }
    }else{
        if(mLeftWriter){
            char* ldata = new char[outstr.size()];
            std::memcpy(ldata, outstr.c_str(), outstr.size());
            mLeftWriter->input(ldata, outstr.size());
        }
    }
    if(mFailedWriter && !failedOut.empty()){
        char* fdata = new char[failedOut.size()];
        std::memcpy(fdata, failedOut.c_str(), failedOut.size());
        mFailedWriter->input(fdata, failedOut.size());
    }
    if(!mOptions->split.enabled){
        mOutputMtx.unlock();
    }
    if(mOptions->split.byFileLines){
        config->markProcessed(readPassed);
    }else{
        config->markProcessed(pack->count);
    }
    delete[] pack->data;
    delete pack;
}

void SingleEndProcessor::writeTask(WriterThread* config){
    while(true){
        if(config->isCompleted()){
            config->output();
            break;
        }
        config->output();
    }
    util::loginfo(config->getFilename() + " writer finished", mOptions->logmtx);
}
