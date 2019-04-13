#include "filter.h"

int Filter::passFilter(Read* r){
    if(r == NULL || r->length() == 0){
        return COMMONCONST::FAIL_LENGTH;
    }

    int rlen = r->length();
    int lowQualNum = 0;
    int nBaseNum = 0;
    int totalQual = 0;
    const char* seq = r->seq.seqStr.c_str();
    const char* qual = r->quality.c_str();

    if(mOptions->qualFilter.enabled || mOptions->lengthFilter.enabled){
        for(int i = 0; i < rlen; ++i){
            totalQual += qual[i] - 33;
            if(seq[i] == 'N'){
                ++nBaseNum;
            }
            if(qual[i] < mOptions->qualFilter.lowQualityLimit){
                ++lowQualNum;
            }
        }
    }
    if(mOptions->qualFilter.enabled){
        if(lowQualNum > mOptions->qualFilter.lowQualityBaseLimit){
            return COMMONCONST::FAIL_QUALITY;
        }else if(mOptions->qualFilter.averageQualityLimit > 0 && mOptions->qualFilter.averageQualityLimit > double(totalQual) / rlen){
            return COMMONCONST::FAIL_QUALITY;
        }
    }

    if(mOptions->qualFilter.enabled && nBaseNum > mOptions->qualFilter.nBaseLimit){
        return COMMONCONST::FAIL_N_BASE;
    }

    if(mOptions->lengthFilter.enabled){
        if(rlen < mOptions->lengthFilter.minReadLength){
            return COMMONCONST::FAIL_LENGTH;
        }
        if(mOptions->lengthFilter.maxReadLength > 0 && rlen > mOptions->lengthFilter.maxReadLength){
            return COMMONCONST::FAIL_TOO_LONG;
        }
    }

    if(mOptions->complexityFilter.enabled && !passLowComplexityFliter(r)){
        return COMMONCONST::FAIL_COMPLEXITY;
    }

    return COMMONCONST::PASS_FILTER;
}

bool Filter::passLowComplexityFliter(Read* r){
    int diff = 0;
    int rlen = r->length();
    if(rlen <= 1){
        return false;
    }
    const char* seq = r->seq.seqStr.c_str();
    for(int i = 0; i < rlen - 1; ++i){
        if(seq[i] != seq[i+1]){
            ++diff;
        }
    }
    return (double)diff/(rlen - 1) >= mOptions->complexityFilter.threshold;
}

Read* Filter::trimAndCut(Read* r, int forceFrontCut, int forceTailCut){
    // do not need quality cutting
    if(forceFrontCut == 0 && forceTailCut == 0 && !mOptions->qualitycut.enableFront && !mOptions->qualitycut.enableRright && !mOptions->qualitycut.enableTail){
        return r;
    }

    int rlen = r->length() - forceFrontCut - forceTailCut;
    if(rlen < 0){
        return NULL;
    }

    if(forceFrontCut == 0 && !mOptions->qualitycut.enableFront && !mOptions->qualitycut.enableRright && !mOptions->qualitycut.enableTail){
        r->resize(rlen);
        return r;
    }else if(!mOptions->qualitycut.enableFront && !mOptions->qualitycut.enableRright && !mOptions->qualitycut.enableTail){
        r->seq.seqStr = r->seq.seqStr.substr(forceFrontCut, rlen);
        r->quality = r->quality.substr(forceFrontCut, rlen);
        return r;
    }

    // need quality cutting
    int l = r->length();
    const char* qual = r->quality.c_str();
    const char* seq = r->seq.seqStr.c_str();
    // quality cutting forward by sliding window
    if(mOptions->qualitycut.enableFront){
        int w = mOptions->qualitycut.windowSizeFront;
        int s = forceFrontCut;
        if(l - forceFrontCut - forceTailCut - w <= 0){
            return NULL;
        }
        int totalQual = 0;
        for(int i = 0; i < w - 1; ++i){
            totalQual += qual[s + i];
        }
        for(s = forceFrontCut; s + w < l - forceTailCut; ++s){
            totalQual += qual[s + w -1];
            if(s > forceFrontCut){
                totalQual -= qual[s - 1];
            }
            if((double)totalQual/(double)w >= 33 + mOptions->qualitycut.qualityFront){
                break;
            }
        }
        if(s > 0){
            s = s + w -1;
        }
        // forceFrontCut is relocated one base before the last sliding window's last base and skip any N afterwards aswell
        while(s < l && seq[s] == 'N'){
            ++s;
        }
        forceFrontCut = s;
        rlen = l - forceFrontCut - forceTailCut;
    }

    // quality cutting from right by sliding window
    if(mOptions->qualitycut.enableRright){
        int w = mOptions->qualitycut.windowSizeRight;
        int s = forceFrontCut;
        if(l - forceFrontCut - forceTailCut - w <= 0){
            return NULL;
        }
        int totalQual = 0;
        bool foundLowQualWindow = false;
        for(int i = 0; i < w - 1; ++i){
            totalQual += qual[s + i];
        }
        for(s = forceFrontCut; s + w < l - forceTailCut; ++s){
            totalQual += qual[s + w -1];
            if(s > forceFrontCut){
                totalQual -= qual[s - 1];
            }
            if((double)totalQual/(double)w < 33 + mOptions->qualitycut.qualityRight){
                foundLowQualWindow = true;
                break;
            }
        }
        if(foundLowQualWindow){
            while(s < l - 1 && qual[s] >= 33 + mOptions->qualitycut.qualityRight){
                ++s;
            }
            rlen = s - forceFrontCut;
        }
    }

    // quality cutting backward by sliding window
    if(!mOptions->qualitycut.enableRright && mOptions->qualitycut.enableTail){
        int w = mOptions->qualitycut.windowSizeTail;
        if(l - forceFrontCut - forceTailCut - w <= 0){
            return NULL;
        }
        int totalQual = 0;
        int t = l - forceTailCut - 1;
        for(int i = 0; i < w - 1; ++i){
            totalQual += qual[t - i];
        }
        for(t = l - forceTailCut - 1; t - w >= forceFrontCut; --t){
            totalQual += qual[t - w + 1];
            if(t < l - forceTailCut - 1){
                totalQual -= qual[t + 1];
            }
            if((double)totalQual/(double)w >= 33 +  mOptions->qualitycut.qualityTail){
                break;
            }
        }
        if(t < l - 1){
            t = t - w + 1;
        }
        while(t >= 0 && seq[t] == 'N'){
            --t;
        }
        rlen = t - forceFrontCut + 1;
    }

    if(rlen <= 0 || forceFrontCut >= l - 1){
        return NULL;
    }
    r->seq.seqStr = r->seq.seqStr.substr(forceFrontCut, rlen);
    r->quality = r->quality.substr(forceFrontCut, rlen);
    return r;
}

bool Filter::match(const std::vector<std::string>& list, const std::string& target, int threshold){
    int tlen = target.length();
    int slen = 0;
    int diff = 0;
    for(size_t i = 0; i < list.size(); ++i){
        diff = 0;
        slen = list[i].length();
        for(int s = 0; s < slen && s < tlen; ++s){
            if(list[i][s] != target[s]){
                ++diff;
                if(diff > threshold){
                    break;
                }
            }
        }
        if(diff <= threshold){
            return true;
        }
    }
    return false;
}

bool Filter::filterByIndex(Read* r){
    if(mOptions->indexFilter.enabled){
        if(match(mOptions->indexFilter.blacklist1, r->firstIndex(), mOptions->indexFilter.threshold)){
            return true;
        }
    }
    return false;
}

bool Filter::filterByIndex(Read* r1, Read* r2){
    if(mOptions->indexFilter.enabled){
        if(match(mOptions->indexFilter.blacklist1, r1->firstIndex(), mOptions->indexFilter.threshold)){
            return true;
        }
        if(match(mOptions->indexFilter.blacklist2, r2->firstIndex(), mOptions->indexFilter.threshold)){
            return true;
        }
    }
    return false;
}
