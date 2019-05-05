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
        usleep(1000);
    }
    mRepo.packBuffer[mRepo.writePos] = pack;
    ++mRepo.writePos;
}

void SingleEndProcessor::producerTask(){
    util::loginfo("loading data started", mOptions->logmtx);
    size_t lastReported = 0; // total number of reads have been loaded into memory
    int slept = 0;  // #sleep happened
    size_t readNum = 0;  // total number of reads have been loaded into memory and put into mReop
    Read** data = new Read*[mOptions->bufSize.maxReadsInPack];
    std::memset(data, 0, sizeof(Read*) * mOptions->bufSize.maxReadsInPack);
    FqReader reader(mOptions->in1, true, mOptions->phred64);
    size_t count = 0; // number of reads have been loaded into memory and put into data but not mReop 
    bool needToBreak = false;
    while(true){
        Read* read = reader.read();
        if(!read || needToBreak){
            // last read have got or the first N reads needed have got
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
        // if only need to process the first N reads
        if(mOptions->readsToProcess > 0 && count + readNum >= mOptions->readsToProcess){
            needToBreak = true;
        }
        if(count + readNum >= lastReported + 1000000){
            lastReported = count + readNum;
            std::string msg = "loaded " + std::to_string(lastReported/1000000) + "M reads";
            util::loginfo(msg, mOptions->logmtx);
        }
        // the pack is full or first N reads needed have got
        if(count == mOptions->bufSize.maxReadsInPack || needToBreak){
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            //re-initialize data for next pack
            data = new Read*[mOptions->bufSize.maxReadsInPack];
            std::memset(data, 0, sizeof(Read*) * mOptions->bufSize.maxReadsInPack);
            // if the consumer is far behind this producer, sleep and wait to limit memory usage
            while(mRepo.writePos - mRepo.readPos > mOptions->bufSize.maxPacksInMemory){
                ++slept;
                usleep(100);
            }
            readNum += count;
            // if the writer threads are far behind this producer, sleep and wait
            if(readNum % (mOptions->bufSize.maxReadsInPack * mOptions->bufSize.maxPacksInMemory) == 0 && mLeftWriter){
                while(mLeftWriter->bufferLength() > mOptions->bufSize.maxPacksInMemory){
                    ++slept;
                    usleep(100);
                }
            }
            // rest count to 0
            count = 0;
        }
    }
    mProduceFinished = true;
    util::loginfo("all reads loaded, start to monitor thread status", mOptions->logmtx);
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
            usleep(1000);
        }
        if(mProduceFinished && mRepo.writePos == mRepo.readPos){
            ++mFinishedThreads;
            std::string msg = "thread " + std::to_string(config->getThreadId() + 1) + " data processing finished";
            util::loginfo(msg, mOptions->logmtx);
            break;
        }
        std::string msg = "thread " + std::to_string(config->getThreadId() + 1) + " is processing the " +
                           std::to_string(mRepo.readPos) + "/" + std::to_string(mRepo.writePos) + " pack";
        util::loginfo(msg, mOptions->logmtx);
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
    std::string msg = "thread " + std::to_string(config->getThreadId() + 1) + " finished";
    util::loginfo(msg, mOptions->logmtx);
}

void SingleEndProcessor::consumePack(ThreadConfig* config){
    ReadPack* data;
    mInputMtx.lock();
    while(mRepo.writePos <= mRepo.readPos){
        usleep(1000);
        if(mProduceFinished){
            mInputMtx.unlock();
            return;
        }
    }
    data = mRepo.packBuffer[mRepo.readPos];
    ++mRepo.readPos;
    if(mRepo.readPos == mOptions->bufSize.maxPacksInReadPackRepo){
        processSingleEnd(data, config);
        mRepo.readPos = mRepo.readPos % mOptions->bufSize.maxPacksInReadPackRepo;
        mRepo.writePos = mRepo.writePos % mOptions->bufSize.maxPacksInReadPackRepo;
        mInputMtx.unlock();
    }else{
        mInputMtx.unlock();
        processSingleEnd(data, config);
    }
}

bool SingleEndProcessor::process(){
    if(!mOptions->split.enabled){
        initOutput();
    }
    initReadPackRepository();
    std::thread producer(std::bind(&SingleEndProcessor::producerTask, this));

    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t = 0; t < mOptions->thread; ++t){
        configs[t] = new ThreadConfig(mOptions, t, false);
        initConfig(configs[t]);
    }

    std::thread** threads = new std::thread*[mOptions->thread];
    for(int t = 0; t < mOptions->thread; ++t){
        threads[t] = new std::thread(std::bind(&SingleEndProcessor::consumerTask, this, configs[t]));
    }

    std::thread* leftWriterThread = NULL;
    if(mLeftWriter){
        leftWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mLeftWriter));
    }
    std::thread* failedWriterThread = NULL;
    if(mFailedWriter){
        failedWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mFailedWriter));
    }
    producer.join();
    for(int t = 0; t < mOptions->thread; ++t){
        threads[t]->join();
    }

    if(!mOptions->split.enabled){
        if(leftWriterThread){
            leftWriterThread->join();
        }
    }
    util::loginfo("start to generate reports\n", mOptions->logmtx);
    // merge stats and read filter results
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
    std::cerr << "Read1 before filtering: " << std::endl;
    std::cerr << finalPreStats << std::endl;
    std::cerr << "Read1 after filtering: " << std::endl;
    std::cerr << finalPostStats << std::endl;
    std::cerr << "Filtering result:" << std::endl;
    std::cerr << finalFilterResult << std::endl;
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
        std::cerr << std::endl;
        std::cerr << "Duplicate rate(may be overestimated since this is SE data): " << dupRate << std::endl;
    }
    
    JsonReporter jr(mOptions);
    jr.setDupHist(dupHist, dupMeanGC, dupRate);
    jr.report(finalFilterResult, finalPreStats, finalPostStats);
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.report(finalFilterResult, finalPreStats, finalPostStats);

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
                PolyX::trimPolyG(r1, mOptions->polyGTrim.minLen);
            }
        }
        // adapter trimming
        if(r1 != NULL && mOptions->adapter.enableTriming && mOptions->adapter.adapterSeqR1Provided){
            AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.inputAdapterSeqR1);
        }
        // polyX trimming
        if(r1 != NULL){
            if(mOptions->polyXTrim.enabled){
                PolyX::trimPolyX(r1, mOptions->polyXTrim.minLen);
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
    std::string msg = config->getFilename() + " writer finished";
    util::loginfo(msg, mOptions->logmtx);
}
