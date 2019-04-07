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
    if(mOptions->out1.empty() || mOptions->out2.empty()){
        return;
    }
    mLeftWriter = new WriterThread(mOptions, mOptions->out1);
    mRightWriter = new WriterThread(mOptions, mOptions->out2);
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
}

void PairEndProcessor::initConfig(ThreadConfig* config){
    if(mOptions->out1.empty()){
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
    std::thread producer(&PairEndProcessor::producerTask, this);
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t = 0; t < mOptions->thread; ++t){
        configs[t] = new ThreadConfig(mOptions, t, true);
        initConfig(configs[t]);
    }
    std::thread** threads = new std::thread*[mOptions->thread];
    for(int t = 0; t < mOptions->thread; ++t){
        threads[t] = new std::thread(&PairEndProcessor::consumePack, this, configs[t]);
    }
    std::thread* leftWriterThread = NULL;
    std::thread* rightWriterThread = NULL;
    if(mLeftWriter){
        leftWriterThread = new std::thread(&PairEndProcessor::writeTask, this, mLeftWriter);
    }
    if(mRightWriter){
        rightWriterThread = new std::thread(&PairEndProcessor::writeTask, this, mRightWriter);
    }
    producer.join();
    for(int t = 0; t < mOptions->thread; ++t){
        threads[t]->join();
    }
    if(leftWriterThread){
        leftWriterThread->join();
    }
    if(rightWriterThread){
        rightWriterThread->join();
    }
    if(mOptions->verbose){
        util::loginfo("start to generate reports\n", mOptions->logmtx);
    }
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

    std::cerr << "Read1 before filtering: \n";
    std::cerr << finalPreStats1 << "\n";
    std::cerr << "Read1 after filtering: \n";
    std::cerr << finalPostStats1 << "\n";
    std::cerr << "Read2 before filtering: \n";
    std::cerr << finalPreStats2 << "\n";
    std::cerr << "Read2 after filtering: \n";
    std::cerr << finalPostStats2 << "\n";
    std::cerr << "Filtering results: \n";
    std::cerr << finalFilterResult << "\n";

    size_t* dupHist = NULL;
    double* dupMeanGC = NULL;
    double dupRate = 0.0;
    if(mOptions->duplicate.enabled){
        dupHist = new size_t[mOptions->duplicate.histSize];
        std::memset(dupHist, 0, sizeof(int) * mOptions->duplicate.histSize);
        dupMeanGC = new double[mOptions->duplicate.histSize];
        std::memset(dupMeanGC, 0, sizeof(double) * mOptions->duplicate.histSize);
        dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
        std::cerr << "\nDuplication rate: " << dupRate * 100.0 << "%\n";
    }

    int peakInsertSize = getPeakInsertSize();
    std::cerr << "Insert size peak (evaluated by pair-end reads: " << peakInsertSize << "\n";
    //clean up
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
    std::string interleaved;
    int readPassed = 0;
    for(int p = 0; p < pack->count; ++p){
        ReadPair* pair = pack->data[p];
        Read* or1 = pair->left;
        Read* or2 = pair->right;
        
        config->getPreStats1()->statRead(or1);
        config->getPreStats2()->statRead(or2);

        if(mDuplicate){
            mDuplicate->statPair(or1, or2);
        }

        if(mOptions->indexFilter.enabled && mFilter->filterByIndex(or1, or2)){
            delete pair;
            continue;
        }

        if(mOptions->umi.enabled){
            mUmiProcessor->process(or1, or2);
        }

        Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1);
        Read* r2 = mFilter->trimAndCut(or2, mOptions->trim.front2, mOptions->trim.tail2);

        if(r1 && r2){
            if(mOptions->polyGTrim.enabled){
                PolyX::trimPolyG(r1, r2, mOptions->polyGTrim.minLen);
            }
        }

        bool insertSizeEvaluated = false;
        if(r1 && r2 && (mOptions->adapter.enableTriming || mOptions->correction.enabled)){
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire);
            if(config->getThreadId() == 0){
                statInsertSize(r1, r2, ov);
                insertSizeEvaluated = true;
            }
            if(mOptions->correction.enabled){
                BaseCorrector::correctByOverlapAnalysis(r1, r2, config->getFilterResult(), ov);
            }
            if(mOptions->adapter.enableTriming){
                bool trimmed = AdapterTrimmer::trimByOverlapAnalysis(r1, r2, config->getFilterResult(), ov);
                if(!trimmed){
                    if(mOptions->adapter.adapterSeqR1Provided){
                        AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.inputAdapterSeqR1, false);
                    }
                    if(mOptions->adapter.adapterSeqR2Provided){
                        AdapterTrimmer::trimBySequence(r2, config->getFilterResult(), mOptions->adapter.inputAdapterSeqR2, true);
                    }
                }
            }

            if(config->getThreadId() == 0 && !insertSizeEvaluated && r1 && r2){
                OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire);
                statInsertSize(r1, r2, ov);
                insertSizeEvaluated = true;
            }

            if(r1 && r2){
                if(mOptions->polyGTrim.enabled){
                    PolyX::trimPolyX(r1,  r2, mOptions->polyXTrim.minLen);
                }
            }

            if(r1 && r2){
                if(mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length()){
                    r1->resize(mOptions->trim.maxLen1);
                }
                if(mOptions->trim.maxLen2 > 0 && mOptions->trim.maxLen2 < r2->length()){
                    r2->resize(mOptions->trim.maxLen2);
                }
            }

            int result1 = mFilter->passFilter(r1);
            int result2 = mFilter->passFilter(r2);
            config->addFilterResult(std::max(result1, result2));

            if(r1 && result1 == COMMONCONST::PASS_FILTER && r2 && result2 == COMMONCONST::PASS_FILTER){
                if(mOptions->outputToSTDOUT){
                    interleaved += r1->toString() + r2->toString();
                }else{
                    outstr1 += r1->toString();
                    outstr2 += r2->toString();
                }
                config->getPostStats1()->statRead(r1);
                config->getPostStats2()->statRead(r2);
                ++readPassed;
            }

            delete pair;
            if(r1 && r1 != or1){
                delete r1;
            }
            if(r2 && r2 != or2){
                delete r2;
            }
        }
    }

    if(!mOptions->split.enabled){
        mOutputMtx.lock();
    }
    if(mOptions->outputToSTDOUT){
        std::fwrite(interleaved.c_str(), 1, interleaved.size(), stdout);
    }else if(mOptions->split.enabled){
        if(!mOptions->out1.empty()){
            config->getWriter1()->writeString(outstr1);
        }
        if(!mOptions->out2.empty()){
            config->getWriter2()->writeString(outstr2);
        }
    }else{
        if(mRightWriter && mLeftWriter){
            char* ldata = new char[outstr1.size()];
            std::memcpy(ldata, outstr1.c_str(), outstr1.size());
            mLeftWriter->input(ldata, outstr1.size());
            char* rdata = new char[outstr2.size()];
            std::memcpy(rdata, outstr2.c_str(), outstr2.size());
            mRightWriter->input(rdata, outstr2.size());
        }
    }
    if(!mOptions->split.enabled){
        mOutputMtx.unlock();
    }
    if(mOptions->split.byFileLines){
        config->markProcessed(readPassed);
    }else{
        config->markProcessed(pack->count);
    }
    delete pack->data;
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
    while(mRepo.writePos >= mOptions->bufSize.maxPacksInMemory){
        sleep(60);
    }
    mRepo.packBuffer[mRepo.writePos] = pack;
    ++mRepo.writePos;
}

void PairEndProcessor::consumePack(ThreadConfig* config){
    mInputMtx.lock();
    while(mRepo.writePos <= mRepo.readPos){
        usleep(1000);
        if(mProduceFinished){
            mInputMtx.unlock();
            return;
        }
    }
    ReadPairPack* data = mRepo.packBuffer[mRepo.readPos];
    ++mRepo.readPos;
    if(mRepo.readPos >= mOptions->bufSize.maxPacksInReadPackRepo){
        mRepo.readPos = mRepo.readPos %  mOptions->bufSize.maxPacksInReadPackRepo;
        mRepo.writePos = mRepo.writePos % mOptions->bufSize.maxPacksInReadPackRepo;
        processPairEnd(data, config);
        mInputMtx.unlock();
    }else{
        mInputMtx.unlock();
        processPairEnd(data, config);
    }
}

void PairEndProcessor::producerTask(){
    if(mOptions->verbose){
        util::loginfo("start to load data", mOptions->logmtx);
    }
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    ReadPair** data = new ReadPair*[mOptions->bufSize.maxReadsInPack];
    std::memset(data, 0, sizeof(ReadPair*) * mOptions->bufSize.maxReadsInPack);
    std::cout << "in1:" << mOptions->in1 << "\n";
    std::cout << "in2:" << mOptions->in2 << "\n";
    FqReaderPair reader(mOptions->in1, mOptions->in2, true, mOptions->phred64, mOptions->interleavedInput);
    int count = 0;
    bool needToBreak = false;
    while(true){
        ReadPair* readPair = reader.read();
        if(!readPair || needToBreak){
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
        if(mOptions->readsToProcess > 0 && count + readNum >= mOptions->readsToProcess){
            needToBreak = true;
        }
        if(mOptions->verbose && count + readNum >= lastReported + 1000000){
            lastReported = count + readNum;
            std::string msg = "loaded " + std::to_string(lastReported/1000000) + "M";
            util::loginfo(msg, mOptions->logmtx);
        }
        if(count == mOptions->bufSize.maxPacksInReadPackRepo || needToBreak){
            ReadPairPack* pack = new ReadPairPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            data = new ReadPair*[mOptions->bufSize.maxPacksInReadPackRepo];
            std::memset(data, 0, sizeof(ReadPair*) * mOptions->bufSize.maxPacksInReadPackRepo);
            while(mRepo.writePos - mRepo.readPos > mOptions->bufSize.maxPacksInMemory){
                ++slept;
                usleep(1000);
            }
            readNum += count;
            if(readNum % (mOptions->bufSize.maxReadsInPack * mOptions->bufSize.maxPacksInMemory) == 0 && mLeftWriter){
                while( (mLeftWriter && mLeftWriter->bufferLength() > mOptions->bufSize.maxPacksInMemory) ||
                       (mRightWriter && mRightWriter->bufferLength() > mOptions->bufSize.maxPacksInMemory)){
                    ++slept;
                    usleep(1000);
                }
            }
            count = 0;
        }
        mProduceFinished = true;
    }
    if(mOptions->verbose){
        util::loginfo("all reads loaded, start to monitor thread status", mOptions->logmtx);
    }
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
            usleep(1000);
        }
        if(mProduceFinished && mRepo.writePos == mRepo.readPos){
            ++mFinishedThreads;
            if(mOptions->verbose){
                std::string msg = "thread " + std::to_string(config->getThreadId() + 1) + " data processing completed";
                util::loginfo(msg, mOptions->logmtx);
            }
            break;
        }
        if(mProduceFinished){
            if(mOptions->verbose){
                std::string msg = "thread " + std::to_string(config->getThreadId() + 1) + " is processing the " +
                                  std::to_string(mRepo.readPos) + " / " + std::to_string(mRepo.writePos) + " pack";
                util::loginfo(msg, mOptions->logmtx);
            }
            consumePack(config);
        }
    }
    if(mFinishedThreads == mOptions->thread){
        if(mLeftWriter){
            mLeftWriter->setInputCompleted();
        }
        if(mRightWriter){
            mRightWriter->setInputCompleted();
        }
    }

    if(mOptions->verbose){
        std::string msg = "thread " + std::to_string(config->getThreadId() + 1) + " finished";
        util::loginfo(msg, mOptions->logmtx);
    }
}

void PairEndProcessor::writeTask(WriterThread* config){
    while(true){
        if(config->isCompleted()){
            config->output();
            break;
        }
        config->output();
    }
    if(mOptions->verbose){
        std::string msg = config->getFilename() + " writer finished";
        util::loginfo(msg, mOptions->logmtx);
    }
}
