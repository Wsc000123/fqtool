#include "polyx.h"

PolyX::PolyX(){
}

PolyX::~PolyX(){
}

void PolyX::trimPolyG(Read* r1, Read* r2, int compareReq, int maxMismatch, int allowedOneMismatchForEach, FilterResult* fr){
    trimPolyG(r1, compareReq, maxMismatch, allowedOneMismatchForEach, fr);
    trimPolyG(r2, compareReq, maxMismatch, allowedOneMismatchForEach, fr);
}

void PolyX::trimPolyG(Read* r, int compareReq, int maxMismatch, int allowedOneMismatchForEach, FilterResult* fr){
    const char* data = r->seq.seqStr.c_str();
    int rlen = r->length();
    int mismatch = 0;
    int i = 0;
    int firstGpos = rlen - 1;
    for(i = 0; i < rlen; ++i){
        if(data[rlen - i - 1] != 'G'){
            ++mismatch;
        }else{
            firstGpos = rlen - i - 1;
        }
        int allowedMismatch = std::min(maxMismatch, std::max(1, (i+1)/allowedOneMismatchForEach));
        if(mismatch > allowedMismatch){
            break;
        }
    }
    if(i + 1 >= compareReq){
        r->resize(firstGpos);
        if(fr){
            fr->addPolyXTrimmed(3, rlen - firstGpos);
        }
    }

}

void PolyX::trimPolyX(Read* r1, Read* r2, std::string trimChr, int compareReq, int maxMismatch, int allowedOneMismatchForEach, FilterResult* fr){
    trimPolyX(r1, trimChr, compareReq, maxMismatch, allowedOneMismatchForEach, fr);
    trimPolyX(r2, trimChr, compareReq, maxMismatch, allowedOneMismatchForEach, fr);
}

void PolyX::trimPolyX(Read* r, std::string trimChr, int compareReq, int maxMismatch, int allowedOneMismatchForEach, FilterResult* fr){
    const char* data = r->seq.seqStr.c_str();
    int rlen = r->length();
    int atcgNumbers[5] = {0, 0, 0, 0, 0};
    char atcgBases[5] = {'A', 'T', 'C', 'G', 'N'};
    int pos = 0;
    for(pos = 0; pos < rlen; ++pos){
        switch(data[rlen - 1 -pos]){
            case 'A':
                ++atcgNumbers[0];
                break;
            case 'T':
                ++atcgNumbers[1];
                break;
            case 'C':
                ++atcgNumbers[2];
                break;
            case 'G':
                ++atcgNumbers[3];
                break;
            default:
                ++atcgNumbers[4];
                break;
        }
        int cmp = (pos + 1);
        int allowedMismatch = std::min(maxMismatch, std::max(1, cmp/allowedOneMismatchForEach));
        bool needToBreak = true;
        for(int b = 0; b < 5; ++b){
            if(trimChr.find(atcgBases[b]) != std::string::npos && cmp - atcgNumbers[b] <= allowedMismatch){
                needToBreak = false;
            }
        }
        if(needToBreak){
            break;
        }
    }

    if(pos + 1 >= compareReq){
        int poly = 0;
        int maxCount = -1;
        for(int b = 0; b < 5; ++b){
            if(trimChr.find(atcgBases[b]) != std::string::npos && atcgNumbers[b] > maxCount){
                maxCount = atcgNumbers[b];
                poly = b;
            }
        }
        char polyBase = atcgBases[poly];
        pos = std::min(rlen - 1, pos);
        while(data[rlen - pos - 1] != polyBase && pos > 0){
            --pos;
        }
        std::cerr << r->name << ":\n" << r->seq.seqStr << "\n" << r->seq.seqStr.substr(rlen - pos - 1) << std::endl << std::endl;
        r->resize(rlen - pos - 1);
        if(fr){
            fr->addPolyXTrimmed(poly, pos + 1);
        }
    }
}
