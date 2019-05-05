#include "duplicate.h"

Duplicate::Duplicate(Options* opt){
    mOptions = opt;
    mKeyLenInBase = mOptions->duplicate.keylen;
    mKeyLenInBit = (1UL << (2 * mKeyLenInBase));
    mDups = new uint64_t[mKeyLenInBit];
    std::memset(mDups, 0, sizeof(uint64_t) * mKeyLenInBit);
    mCounts = new uint32_t[mKeyLenInBit];
    std::memset(mCounts, 0, sizeof(uint32_t) * mKeyLenInBit);
    mGC = new uint8_t[mKeyLenInBit];
    std::memset(mGC, 0, sizeof(uint8_t) * mKeyLenInBit);
}

Duplicate::~Duplicate(){
    delete[] mDups;
    delete[] mCounts;
    delete[] mGC;
}

uint64_t Duplicate::seq2int(const char* cstr, int start, int keylen, bool& valid){
    uint64_t ret = 0;
    for(int i = 0; i < keylen; ++i){
        ret <<= 2;
        switch(cstr[start + i]){
            case 'A':
                ret += 0;
                break;
            case 'T':
                ret += 1;
                break;
            case 'C':
                ret += 2;
                break;
            case 'G':
                ret += 3;
                break;
            default:
                valid = false;
                return 0;
        }
    }
    return ret;
}

void Duplicate::addRecord(uint32_t key, uint64_t kmer32, uint8_t gc){
    mAddLock.lock();
    if(mCounts[key] == 0){
        mCounts[key] = 1;
        mDups[key] = kmer32;
        mGC[key] = gc;
    }else{
        if(mDups[key] == kmer32){
            ++mCounts[key];
        }else if(mDups[key] > kmer32){
            mDups[key] = kmer32;
            mCounts[key] = 1;
            mGC[key] = gc;
        }
    }
    mAddLock.unlock();
}

void Duplicate::statRead(Read* r){
    if(r->length() < 32){
        return;
    }

    int start1 = 0;
    int start2 = std::max(0, r->length() - 32 - 5);

    const char* cstr = r->seq.seqStr.c_str();
    bool valid = true;
    uint64_t ret = seq2int(cstr, start1, mKeyLenInBase, valid);
    uint32_t key = (uint32_t)ret;
    if(!valid){
        return;
    }
    uint64_t kmer32 = seq2int(cstr, start2, 32, valid);
    if(!valid){
        return;
    }
    uint8_t gc = 0;
    if(mCounts[key] == 0){
        for(int i = 0; i < r->length(); ++i){
            if(cstr[i] == 'C' || cstr[i] == 'G'){
                ++gc;
            }
        }
    }
    gc = std::round(255.0 * (double) gc / (double) r->length());
    addRecord(key, kmer32, gc);
}

void Duplicate::statPair(Read* r1, Read* r2){
    if(r1->length() < 32 || r2->length() < 32){
        return;
    }
    const char* cstr1 = r1->seq.seqStr.c_str();
    const char* cstr2 = r2->seq.seqStr.c_str();
    bool valid = true;

    uint64_t ret = seq2int(cstr1, 0, mKeyLenInBase, valid);
    uint32_t key = (uint32_t)ret;
    if(!valid){
        return;
    }
    uint64_t kmer32 = seq2int(cstr2, 0, 32, valid);
    if(!valid){
        return;
    }

    uint8_t gc = 0;
    if(mCounts[key] == 0){
        for(int i = 0; i < r1->length(); ++i){
            if(cstr1[i] == 'C' || cstr1[i] == 'G'){
                ++gc;
            }
        }
        for(int i = 0; i < r2->length(); ++i){
            if(cstr2[i] == 'C' || cstr2[i] == 'G'){
                ++gc;
            }
        }
    }

    gc = std::round(255.0 * (double) gc / (double)(r1->length() + r2->length()));
    addRecord(key, kmer32, gc);
}

double Duplicate::statAll(size_t* hist, double* meanGC, size_t histSize){
    size_t totalNum = 0;
    size_t dupNum = 0;
    int* gcStatNum = new int[histSize];
    std::memset(gcStatNum, 0, sizeof(int) * histSize);
    for(uint64_t key = 0; key < mKeyLenInBit; ++key){
        uint32_t count = mCounts[key];
        uint8_t gc = mGC[key];
        if(count > 0){
            totalNum += count;
            dupNum += count - 1;
            if(count > histSize){
                ++hist[histSize - 1];
                meanGC[histSize - 1] += gc;
                ++gcStatNum[histSize - 1];
            }else{
                ++hist[count];
                meanGC[count] += gc;
                ++gcStatNum[count];
            }
        }
    }
    
    for(size_t i = 0; i < histSize; ++i){
        if(gcStatNum[i] > 0){
            meanGC[i] = meanGC[i] / 255.0 / gcStatNum[i];
        }
    }

    delete[] gcStatNum;
    if(totalNum == 0){
        return 0.0;
    }else{
        return (double)dupNum / (double)totalNum;
    }
}
