#include "adaptertrimmer.h"

AdapterTrimmer::AdapterTrimmer(){
}

AdapterTrimmer::~AdapterTrimmer(){
}

bool AdapterTrimmer::trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr){
    OverlapResult ov = OverlapAnalysis::analyze(r1, r2);
    return trimByOverlapAnalysis(r1, r2, fr, ov);
}

bool AdapterTrimmer::trimByOverlapAnalysis(Read* r1, Read* r2, FilterResult* fr, OverlapResult& ov){
    int ol = ov.overlapLen;
    if(ov.diff <= 5 && ov.overlapped && ov.offset < 0 && ol > r1->length() / 3){
        std::string adapter1 = r1->seq.seqStr.substr(ol, r1->length() - ol);
        std::string adapter2 = r2->seq.seqStr.substr(ol, r2->length() - ol);
        r1->seq.seqStr = r1->seq.seqStr.substr(0, ol);
        r1->quality = r1->quality.substr(0, ol);
        r2->seq.seqStr = r2->seq.seqStr.substr(0, ol);
        r2->quality = r2->quality.substr(0, ol);
        fr->addAdapterTrimmed(adapter1, adapter2);
        return true;
    }
    return false;
}

bool AdapterTrimmer::trimBySequence(Read* r, FilterResult* fr, std::string& adapterSeq, bool isR2){
    const int matchRequired = 4;
    const int allowOneMismatchForEach = 8;
    int rlen = r->length();
    int alen = adapterSeq.length();
    const char* rdata = r->seq.seqStr.c_str();
    const char* adata = adapterSeq.c_str();

    if(alen < matchRequired){
        return false;
    }

    int pos = 0;
    bool found = false;
    int start = 0;
    // skip first few adapter sequences as they might be polyA 
    if(alen >= 16){
        start = -4;
    }else if(alen >= 12){
        start -= 3;
    }else if(alen >= 8){
        start = -2;
    }
    for(pos = start; pos < rlen - matchRequired; ++pos){
        int cmplen = std::min(rlen - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int mismatch = 0;
        bool matched = true;
        for(int i = std::max(0, -pos); i < cmplen; ++i){
            if(adata[i] != rdata[i + pos]){
                ++mismatch;
                if(mismatch > allowedMismatch){
                    matched = false;
                    break;
                }
            }
        }
        if(matched){
            found = true;
            break;
        }
    }
    if(found){
        if(pos < 0){
            std::string adapter = adapterSeq.substr(-pos, alen + pos);
            r->seq.seqStr.resize(0);
            r->quality.resize(0);
            if(fr){
                fr->addAdapterTrimmed(adapter, isR2);
            }
        }else{
            std::string adapter = r->seq.seqStr.substr(pos, rlen - pos);
            r->seq.seqStr = r->seq.seqStr.substr(0, pos);
            r->quality = r->quality.substr(0, pos);
            if(fr){
                fr->addAdapterTrimmed(adapter, isR2);
            }
        }
        return true;
    }
    return false;
}
