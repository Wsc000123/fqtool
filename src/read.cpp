#include "read.h"

Read* ReadPair::merge(){
    Read* rcRight = right->reverseComplement();
    int lenl = left->length();
    int lenr = rcRight->length();
    // use pointer to accelerate operations
    const char* pseql = left->seq.seqStr.c_str();
    const char* pseqr = rcRight->seq.seqStr.c_str();
    const char* pquall = left->quality.c_str();
    const char* pqualr = rcRight->quality.c_str();

    // minimum overlap needed to merge a pair of reads
    const int MIN_OVERLAP = 30;
    bool overlapped = false;
    int olen = MIN_OVERLAP;
    // diff count for a pair of different bases at the same position in an overlap
    int diff = 0;
    // diff count for a pair of low/high quality at the same position in an overlap
    int lowQualDiff = 0;

    while(olen <= std::min(lenl, lenr)){
        diff = 0;
        lowQualDiff = 0;
        bool validOVerlap = true;
        int offset = lenl - olen;
        for(int i = 0; i < olen; ++i){
            if(pseql[offset + i] != pseqr[i]){
                ++diff;
                // lowQualDiff will increase if one base quality >= 30('?') and the other base quality <= 15('0')
                if((pquall[offset + i] >= '?' && pqualr[i] <= '0') || (pquall[offset + i] <= '0' && pqualr[i] >= '?')){
                    ++lowQualDiff;
                }
                // high quality base diff or more than 3 low quality base diff will invalid overlap not shorter than this one
                if(diff > lowQualDiff || lowQualDiff >= 3){
                    validOVerlap = false;
                    break;
                }
            }
        }
        // if found valid overlap, then extend the overlap by one base, else break
        if(validOVerlap){
            overlapped = true;
            break;
        }
        ++olen;
    }

    // if overlapped, merge the reads
    if(overlapped){
        int offset = lenl - olen;
        std::stringstream ss;
        ss << left->name << " merged offset: " << offset << " overlap: " << olen << " diff: " << diff;
        std::string mergedName = ss.str();
        std::string mergedSeq = left->seq.seqStr.substr(0, offset) + rcRight->seq.seqStr;
        std::string mergedQual = left->quality.substr(0, offset) + rcRight->quality;
        // quality adjustion and base calling correction for low quality diff bases
        for(int i = 0; i < olen; ++i){
            // if lowQualDiff happens, keep the base with high quality
            // else use the sum of quality score as the quality 
            if(pseql[offset + i] != pseqr[i]){
                if(pquall[offset + i] >= '?' && pqualr[i] <= '0'){
                    mergedSeq[offset + i] = pseql[offset + i];
                    mergedQual[offset + i] = pquall[offset + i];
                }else{
                    mergedSeq[offset + i] = pseqr[i];
                    mergedQual[offset + i] = pqualr[i];
                }
            }else{
                mergedQual[offset + i] = pquall[offset + i] + pqualr[i] - 33;
            }
        }
        delete rcRight;
        return new Read(mergedName, mergedSeq, "+", mergedQual);
    }
    delete rcRight;
    return NULL;
}
