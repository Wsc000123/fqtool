#include "peprocessor.h"

PairEndProcessor::PairEndProcessor(Options* opt){
    mOptions = opt;
    mProduceFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
    mOutStream1 = NULL;
    mOutStream2 = NULL;
    mZipFile1 = NULL;
    mZipFile2 = NULL;
    mUmiProcessor = new UmiProcessor(mOptions);

    int insertSizeBufLen = mOptions->insertSizeMax + 1;
    mInsertSizeHist = new long[insertSizeBufLen];
    std::memset(mInsertSizeHist, 0, sizeof(long) * insertSizeBufLen);
    mLeftWriter = NULL;
    mRightWriter = NULL;
    mUnPairedLeftWriter = NULL;
    mUnPairedRightWriter = NULL;
    mMergedWriter = NULL;
    mFailedWriter = NULL;

    mDuplicate = NULL;
    if(mOptions->duplicate.enabled){
        mDuplicate = new Duplicate(mOptions);
    }
}

PairEndProcessor::~PairEndProcessor(){
    delete mInsertSizeHist;
    delete mUmiProcessor;
    if(mDuplicate){
        delete mDuplicate;
        mDuplicate = NULL;
    }
}

void PairEndProcessor::initOutput(){
    if(!mOptions->unpaired1.empty()){
        mUnPairedLeftWriter = new WriterThread(mOptions, mOptions->unpaired1);
    }
    if(!mOptions->unpaired2.empty() && mOptions->unpaired2 != mOptions->unpaired1){
        mUnPairedRightWriter = new WriterThread(mOptions, mOptions->unpaired2);
    }
    if(mOptions->mergePE.enabled){
        if(!mOptions->mergePE.out.empty()){
            mMergedWriter = new WriterThread(mOptions, mOptions->mergePE.out);
        }
    }
    if(!mOptions->failedOut.empty()){
        mFailedWriter = new WriterThread(mOptions, mOptions->failedOut);
    }
    if(mOptions->out1.empty()){
        return;
    }
    mLeftWriter = new WriterThread(mOptions, mOptions->out1);
    if(!mOptions->out2.empty()){
        mRightWriter = new WriterThread(mOptions, mOptions->out2);
    }
}

void PairEndProcessor::closeOutput(){
    if(mLeftWriter){
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
    if(mRightWriter){
        delete mRightWriter;
        mRightWriter = NULL;
    }
    if(mMergedWriter){
        delete mMergedWriter;
        mMergedWriter = NULL;
    }
    if(mFailedWriter){
        delete mFailedWriter;
        mFailedWriter = NULL;
    }
    if(mUnPairedLeftWriter){
        delete mUnPairedLeftWriter;
        mUnPairedLeftWriter = NULL;
    }
    if(mUnPairedRightWriter){
        delete mUnPairedRightWriter;
        mUnPairedRightWriter = NULL;
    }
}

void PairEndProcessor::initConfig(ThreadConfig* config){
    if(mOptions->out1.empty() || mOptions->out2.empty()){
        return;
    }
    if(mOptions->split.enabled){
        config->initWriterForSplit();
    }
}

bool PairEndProcessor::process(){
    if(!mOptions->split.enabled){
        initOutput();
    }
    initReadPairPackRepository();
    util::loginfo("read pack repo initialized", mOptions->logmtx);
    std::thread producer(&PairEndProcessor::producerTask, this);
    util::loginfo("producer thread started", mOptions->logmtx);
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t = 0; t < mOptions->thread; ++t){
        configs[t] = new ThreadConfig(mOptions, t, true);
        initConfig(configs[t]);
    }
    std::thread** threads = new std::thread*[mOptions->thread];
    for(int t = 0; t < mOptions->thread; ++t){
        threads[t] = new std::thread(&PairEndProcessor::consumerTask, this, configs[t]);
    }
    util::loginfo(std::to_string(mOptions->thread) + " working threads started", mOptions->logmtx);
    std::thread* leftWriterThread = NULL;
    std::thread* rightWriterThread = NULL;
    std::thread* unpairedLeftWriterThread = NULL;
    std::thread* unpairedRightWriterThread = NULL;
    std::thread* mergedWriterThread = NULL;
    std::thread* failedWriterThread = NULL;
    if(mLeftWriter){
        leftWriterThread = new std::thread(&PairEndProcessor::writeTask, this, mLeftWriter);
         util::loginfo("read1 writer thread started", mOptions->logmtx);
    }
    if(mRightWriter){
        rightWriterThread = new std::thread(&PairEndProcessor::writeTask, this, mRightWriter);
        util::loginfo("read2 writer thread started", mOptions->logmtx);
    }
    if(mUnPairedLeftWriter){
        unpairedLeftWriterThread = new std::thread(&PairEndProcessor::writeTask, this, mUnPairedLeftWriter);
        util::loginfo("unpaired read1 writer thread started", mOptions->logmtx);
    }
    if(mUnPairedRightWriter){
        unpairedRightWriterThread = new std::thread(&PairEndProcessor::writeTask, this, mUnPairedRightWriter);
        util::loginfo("unpaired read2 writer thread started", mOptions->logmtx);
    }
    if(mFailedWriter){
       failedWriterThread  = new std::thread(&PairEndProcessor::writeTask, this, mFailedWriter);
       util::loginfo("failed reads writer thread started", mOptions->logmtx);
    }
    if(mMergedWriter){
        mergedWriterThread = new std::thread(&PairEndProcessor::writeTask, this, mMergedWriter);
        util::loginfo("mreged reads writer thread started", mOptions->logmtx);
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
        if(rightWriterThread){
            rightWriterThread->join();
            util::loginfo("read2 writer thread finished", mOptions->logmtx);
        }
        if(unpairedLeftWriterThread){
            unpairedLeftWriterThread->join();
            util::loginfo("unpaired read1 writer thread finished", mOptions->logmtx);
        }
        if(unpairedRightWriterThread){
            unpairedRightWriterThread->join();
            util::loginfo("unpaired read2 writer thread finished", mOptions->logmtx);
        }
        if(mergedWriterThread){
            mergedWriterThread->join();
            util::loginfo("mreged reads writer thread finished", mOptions->logmtx);
        }
        if(failedWriterThread){
            failedWriterThread->join();
            util::loginfo("failed reads writer thread finished", mOptions->logmtx);
        }
    }
    util::loginfo("start generating reports", mOptions->logmtx);
    std::vector<Stats*> preStats1;
    std::vector<Stats*> preStats2;
    std::vector<Stats*> postStats1;
    std::vector<Stats*> postStats2;
    std::vector<FilterResult*> filterResults;
    for(int t = 0; t < mOptions->thread; ++t){
        preStats1.push_back(configs[t]->getPreStats1());
        preStats2.push_back(configs[t]->getPreStats2());
        postStats1.push_back(configs[t]->getPostStats1());
        postStats2.push_back(configs[t]->getPostStats2());
        filterResults.push_back(configs[t]->getFilterResult());
    }
    Stats* finalPreStats1 = Stats::merge(preStats1);
    Stats* finalPreStats2 = Stats::merge(preStats2);
    Stats* finalPostStats1 = Stats::merge(postStats1);
    Stats* finalPostStats2 = Stats::merge(postStats2);
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);
    // duplication analysis
    size_t* dupHist = NULL;
    double* dupMeanGC = NULL;
    double dupRate = 0.0;
    if(mOptions->duplicate.enabled){
        dupHist = new size_t[mOptions->duplicate.histSize];
        std::memset(dupHist, 0, sizeof(size_t) * mOptions->duplicate.histSize);
        dupMeanGC = new double[mOptions->duplicate.histSize];
        std::memset(dupMeanGC, 0, sizeof(double) * mOptions->duplicate.histSize);
        dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
    }

    int peakInsertSize = getPeakInsertSize();
    JsonReporter jr(mOptions);
    jr.setDupHist(dupHist, dupMeanGC, dupRate);
    jr.setInsertHist(mInsertSizeHist, peakInsertSize);
    jr.report(finalFilterResult,finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.report(finalFilterResult,finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);
    util::loginfo("finish generating reports", mOptions->logmtx);
    // clean up
    for(int t = 0; t < mOptions->thread; ++t){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }
    delete finalPreStats1;
    delete finalPreStats2;
    delete finalPostStats1;
    delete finalPostStats2;
    delete finalFilterResult;
    if(mOptions->duplicate.enabled){
        delete[] dupHist;
        delete[] dupMeanGC;
    }
    delete[] threads;
    delete[] configs;
    if(leftWriterThread){
        delete leftWriterThread;
    }
    if(rightWriterThread){
        delete rightWriterThread;
    }
    if(!mOptions->split.enabled){
        closeOutput();
    }
    return true;
}

int PairEndProcessor::getPeakInsertSize(){
    int peak = 0;
    long maxCount = -1;
    for(int i = 0; i < mOptions->insertSizeMax; ++i){
        if(mInsertSizeHist[i] > maxCount){
            peak = i;
            maxCount = mInsertSizeHist[i];
        }
    }
    return peak;
}

bool PairEndProcessor::processPairEnd(ReadPairPack* pack, ThreadConfig* config){
    std::string outstr1;
    std::string outstr2;
    std::string failedOut;
    std::string unpairedOut1;
    std::string unpairedOut2;
    std::string mergedOutput; 
    std::string singleOutput;
    int readPassed = 0;
    int mergedCount = 0;
    for(int p = 0; p < pack->count; ++p){
        ReadPair* pair = pack->data[p];
        Read* or1 = pair->left;
        Read* or2 = pair->right;
        // do preprocess statistics
        config->getPreStats1()->statRead(or1);
        config->getPreStats2()->statRead(or2);
        // do duplicate analysis if enabled
        if(mOptions->duplicate.enabled){
            mDuplicate->statPair(or1, or2);
        }
        // filter by index if enabled
        if(mOptions->indexFilter.enabled && mFilter->filterByIndex(or1, or2)){
            delete pair;
            continue;
        }
        // process umi if enabled
        if(mOptions->umi.enabled){
            mUmiProcessor->process(or1, or2);
        }
        // trim and cut by length and quality if enabled
        Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1);
        Read* r2 = mFilter->trimAndCut(or2, mOptions->trim.front2, mOptions->trim.tail2);
        // trim polyG if enabled
        if(r1 && r2){
            if(mOptions->polyGTrim.enabled){
                PolyX::trimPolyG(r1, r2, mOptions->polyGTrim.minLen);
            }
        }
        // do insertsize statistics only in thread 0, do adapter trimming and base correction if enabled
        bool insertSizeEvaluated = false;
        if(r1 && r2 && (mOptions->adapter.enableTriming || mOptions->correction.enabled)){
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire);
            // first stat insertsize
            if(config->getThreadId() == 0){
                statInsertSize(r1, r2, ov);
                insertSizeEvaluated = true;
            }
            // second do base correction
            if(mOptions->correction.enabled){
                BaseCorrector::correctByOverlapAnalysis(r1, r2, config->getFilterResult(), ov);
            }
            // then do adapter trimming
            if(mOptions->adapter.enableTriming){
                // trim by overlap analysis firstly
                bool trimmed = AdapterTrimmer::trimByOverlapAnalysis(r1, r2, config->getFilterResult(), ov);
                // if failed, trim by input adapter if possible
                if(!trimmed){
                    if(mOptions->adapter.adapterSeqR1Provided){
                        AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.inputAdapterSeqR1, false);
                    }
                    if(mOptions->adapter.adapterSeqR2Provided){
                        AdapterTrimmer::trimBySequence(r2, config->getFilterResult(), mOptions->adapter.inputAdapterSeqR2, true);
                    }
                }
            }
        }
        // do insertsize statistics if haven't done in previous step
        if(config->getThreadId() == 0 && !insertSizeEvaluated && r1 && r2){
              OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire);
              statInsertSize(r1, r2, ov);
              insertSizeEvaluated = true;
        }
        // trim polyX if enabled
        if(r1 && r2){
            if(mOptions->polyXTrim.enabled){
                PolyX::trimPolyX(r1, r2, mOptions->polyXTrim.minLen);
            }
        }
        // trim too long read if enabled
        if(r1 && r2){
            if(mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length()){
                r1->resize(mOptions->trim.maxLen1);
            }
            if(mOptions->trim.maxLen2 > 0 && mOptions->trim.maxLen2 < r2->length()){
                r2->resize(mOptions->trim.maxLen2);
            }
        }
        // merge reads if enabled
        Read* merged = NULL;
        bool mergeProcessed = false;
        if(mOptions->mergePE.enabled && r1 && r2){
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire);
            if(ov.overlapped){
                merged = OverlapAnalysis::merge(r1, r2, ov);
                int result = mFilter->passFilter(merged);
                config->addFilterResult(result, 2);
                if(result == COMMONCONST::PASS_FILTER){
                    mergedOutput += merged->toString();
                    config->getPostStats1()->statRead(merged);
                    ++readPassed;
                    ++mergedCount;
                }
                delete merged;
                mergeProcessed = true;
            }else if(!mOptions->mergePE.discardUnmerged){
                int result1 = mFilter->passFilter(r1);
                config->addFilterResult(result1, 1);
                if(result1 == COMMONCONST::PASS_FILTER){
                    mergedOutput += r1->toString();
                    config->getPostStats1()->statRead(r1);
                }
                int result2 = mFilter->passFilter(r2);
                config->addFilterResult(result2, 1);
                if(result2 == COMMONCONST::PASS_FILTER){
                    mergedOutput += r2->toString();
                    config->getPostStats2()->statRead(r2);
                }
                if(result1 == COMMONCONST::PASS_FILTER && result2 == COMMONCONST::PASS_FILTER){
                    ++readPassed;
                }
                mergeProcessed = true;
            }
        }

        if(!mergeProcessed){
            int result1 = mFilter->passFilter(r1);
            int result2 = mFilter->passFilter(r2);
            config->addFilterResult(std::max(result1, result2));

            if(r1 && result1 == COMMONCONST::PASS_FILTER && r2 && result2 == COMMONCONST::PASS_FILTER){
                if(mOptions->outputToSTDOUT){
                    singleOutput += r1->toString() + r2->toString();
                }else{
                    outstr1 += r1->toString();
                    outstr2 += r2->toString();
                }
                if(!mOptions->mergePE.enabled){
                    config->getPostStats1()->statRead(r1);
                    config->getPostStats2()->statRead(r2);
                }
                ++readPassed;
            }else if(r1 && result1 == COMMONCONST::PASS_FILTER){
                if(mUnPairedLeftWriter){
                    unpairedOut1 += r1->toString();
                    if(mFailedWriter){
                        failedOut += or2->toStringWithTag(COMMONCONST::FAILED_TYPES[result2]);
                    }
                }else{
                    if(mFailedWriter){
                        failedOut += or1->toStringWithTag("paired_read_is_failing");
                        failedOut += or2->toStringWithTag(COMMONCONST::FAILED_TYPES[result2]);
                    }
                }
            }else if(r2 && result2 == COMMONCONST::PASS_FILTER){
                if(mUnPairedLeftWriter){
                    unpairedOut2 += r2->toString();
                    if(mFailedWriter){
                        failedOut += or1->toStringWithTag(COMMONCONST::FAILED_TYPES[result2]);
                    }
                }else{
                    if(mFailedWriter){
                        failedOut += or1->toStringWithTag(COMMONCONST::FAILED_TYPES[result1]);
                        failedOut += or2->toStringWithTag("paired_read_is_failing");
                    }
                }
            }
        }

        delete pair;
        if(r1 && r1 != or1){
            delete r1;
        }
        if(r2 && r2 != or2){
            delete r2;
        }
    }
    if(!mOptions->split.enabled){
        mOutputMtx.lock();
    }
    if(mOptions->outputToSTDOUT){
        if(mOptions->mergePE.enabled){
            std::fwrite(mergedOutput.c_str(), 1, mergedOutput.length(), stdout);
        }else{
            std::fwrite(singleOutput.c_str(), 1, singleOutput.size(), stdout);
        }
    }else if(mOptions->split.enabled){
        if(!mOptions->out1.empty()){
            config->getWriter1()->writeString(outstr1);
        }
        if(!mOptions->out2.empty()){
            config->getWriter2()->writeString(outstr2);
        }
    }
    
    if(mMergedWriter && !mergedOutput.empty()){
        char* mdata = new char[mergedOutput.size()];
        std::memcpy(mdata, mergedOutput.c_str(), mergedOutput.size());
        mMergedWriter->input(mdata, mergedOutput.size());
    }

    if(mFailedWriter && !failedOut.empty()){
        char* fdata = new char[failedOut.size()];
        std::memcpy(fdata, failedOut.c_str(), failedOut.size());
        mFailedWriter->input(fdata, failedOut.size());
    }

    if(mRightWriter && mLeftWriter && (!outstr1.empty() || !outstr2.empty())){
        char* ldata = new char[outstr1.size()];
        std::memcpy(ldata, outstr1.c_str(), outstr1.size());
        mLeftWriter->input(ldata, outstr1.size());
        char* rdata = new char[outstr2.size()];
        std::memcpy(rdata, outstr2.c_str(), outstr2.size());
        mRightWriter->input(rdata, outstr2.size());
    }else if(mLeftWriter && !singleOutput.empty()){
        char* ldata = new char[singleOutput.size()];
        std::memcpy(ldata, singleOutput.c_str(), singleOutput.size());
        mLeftWriter->input(ldata, singleOutput.size());
    }

    if(!unpairedOut1.empty() && mUnPairedLeftWriter){
        char* unpairedData1 = new char[unpairedOut1.size()];
        std::memcpy(unpairedData1, unpairedOut1.c_str(), unpairedOut1.size());
        mUnPairedLeftWriter->input(unpairedData1, unpairedOut1.size());
    }
    if(!unpairedOut2.empty() && mUnPairedRightWriter){
        char* unpairedData2 = new char[unpairedOut2.size()];
        std::memcpy(unpairedData2, unpairedOut2.c_str(), unpairedOut2.size());
        mUnPairedRightWriter->input(unpairedData2, unpairedOut2.size());
    }

    if(!mOptions->split.enabled){
        mOutputMtx.unlock();
    }
    if(mOptions->split.byFileLines){
        config->markProcessed(readPassed);
    }else{
        config->markProcessed(pack->count);
    }
    if(mOptions->mergePE.enabled){
        config->addMergedPairs(mergedCount);
    }
    delete[] pack->data;
    delete pack;
    
    return true;
}

void PairEndProcessor::statInsertSize(Read* r1, Read* r2, OverlapResult& ov){
    int isize = mOptions->insertSizeMax;
    if(ov.overlapped){
        if(ov.offset > 0){
            isize = r1->length() + r2->length() - ov.overlapLen;
        }else{
            isize = ov.overlapLen;
        }
    }
    if(isize > mOptions->insertSizeMax){
        isize = mOptions->insertSizeMax;
    }
    ++mInsertSizeHist[isize];
}

void PairEndProcessor::initReadPairPackRepository(){
    mRepo.packBuffer = new ReadPairPack*[mOptions->bufSize.maxPacksInReadPackRepo];
    std::memset(mRepo.packBuffer, 0, sizeof(ReadPairPack*) * mOptions->bufSize.maxPacksInReadPackRepo);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
}

void PairEndProcessor::destroyReadPairPackRepository(){
    delete[] mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void PairEndProcessor::producePack(ReadPairPack* pack){
    while(mRepo.writePos >= mOptions->bufSize.maxPacksInReadPackRepo){
        usleep(1);
    }
    mRepo.packBuffer[mRepo.writePos] = pack;
    util::loginfo("producer produced pack " + std::to_string(mRepo.writePos), mOptions->logmtx);
    ++mRepo.writePos;
}

void PairEndProcessor::consumePack(ThreadConfig* config){
    mInputMtx.lock();
    while(mRepo.writePos <= mRepo.readPos){
        usleep(1);
        if(mProduceFinished){
            mInputMtx.unlock();
            return;
        }
    }
    size_t packNum = mRepo.readPos;
    ReadPairPack* data = mRepo.packBuffer[packNum];
    ++mRepo.readPos;
    if(mRepo.readPos >= mOptions->bufSize.maxPacksInReadPackRepo){
        mRepo.readPos = mRepo.readPos %  mOptions->bufSize.maxPacksInReadPackRepo;
        mRepo.writePos = mRepo.writePos % mOptions->bufSize.maxPacksInReadPackRepo;
    }
    mInputMtx.unlock();
    util::loginfo("thread " + std::to_string(config->getThreadId()) + " start processing pack " + std::to_string(packNum), mOptions->logmtx);
    processPairEnd(data, config);
    util::loginfo("thread " + std::to_string(config->getThreadId()) + " finish processing pack " + std::to_string(packNum), mOptions->logmtx);
}

void PairEndProcessor::producerTask(){
    util::loginfo("loading data started", mOptions->logmtx);
    size_t readNum = 0;
    ReadPair** data = new ReadPair*[mOptions->bufSize.maxReadsInPack];
    std::memset(data, 0, sizeof(ReadPair*) * mOptions->bufSize.maxReadsInPack);
    FqReaderPair reader(mOptions->in1, mOptions->in2, true, mOptions->phred64, mOptions->interleavedInput);
    size_t count = 0;
    while(true){
        ReadPair* readPair = reader.read();
        ++readNum;
        if(!readPair){
            if(count == 0){
                break;
            }
            ReadPairPack* pack = new ReadPairPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            data = NULL;
            if(readPair){
                delete readPair;
                readPair = NULL;
            }
            break;
        }
        data[count] = readPair;
        ++count;
        if(count == mOptions->bufSize.maxReadsInPack){
            ReadPairPack* pack = new ReadPairPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            data = new ReadPair*[mOptions->bufSize.maxReadsInPack];
            std::memset(data, 0, sizeof(ReadPair*) * mOptions->bufSize.maxReadsInPack);
            while(mRepo.writePos - mRepo.readPos > mOptions->bufSize.maxPacksInMemory){
                usleep(1);
            }
            readNum += count;
            if(readNum % (mOptions->bufSize.maxReadsInPack * mOptions->bufSize.maxPacksInMemory) == 0 && mLeftWriter){
                while( (mLeftWriter && mLeftWriter->bufferLength() > mOptions->bufSize.maxPacksInMemory) ||
                       (mRightWriter && mRightWriter->bufferLength() > mOptions->bufSize.maxPacksInMemory)){
                    usleep(1);
                }
            }
            count = 0;
        }
    }
    mProduceFinished = true;
    util::loginfo("loaded reads: " + std::to_string(readNum - 1), mOptions->logmtx);
}

void PairEndProcessor::consumerTask(ThreadConfig* config){
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
        if(mRightWriter){
            mRightWriter->setInputCompleted();
        }
        if(mUnPairedLeftWriter){
            mUnPairedLeftWriter->setInputCompleted();
        }
        if(mUnPairedRightWriter){
            mUnPairedRightWriter->setInputCompleted();
        }
        if(mMergedWriter){
            mMergedWriter->setInputCompleted();
        }
        if(mFailedWriter){
            mFailedWriter->setInputCompleted();
        }
    }
    util::loginfo("thread " + std::to_string(config->getThreadId()) + " finished", mOptions->logmtx);
}

void PairEndProcessor::writeTask(WriterThread* config){
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
