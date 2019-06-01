#include "fqreader.h"

FqReader::FqReader(const std::string& filename, const bool& hasQuality, const bool& phread64){
    mFileName = filename;
    mGzipFile = NULL;
    mFile = NULL;
    mStdinMode = false;
    mPhread64 = phread64;
    mHasQuality = hasQuality;
    mFqBufSize = (1 << 20);
    mBuf = new char[mFqBufSize];
    mBufDataLen = 0;
    mBufUsedLen = 0;
    mNoLineBreakAtEnd = false;
    init();
}

FqReader::~FqReader(){
    close();
    delete mBuf;
    mBuf = nullptr;
}

bool FqReader::hasNoLineBreakAtEnd(){
    return mNoLineBreakAtEnd;
}

void FqReader::readToBuf(){
    if(isZipped()){
        mBufDataLen = ::gzread(mGzipFile, mBuf, mFqBufSize);
        if(mBufDataLen == -1){
            std::cerr << "Error to read gzip file" << std::endl;
        }
    }else{
        mBufDataLen = std::fread(mBuf, 1, mFqBufSize, mFile);
    }
    mBufUsedLen = 0;
    if(mBufDataLen < mFqBufSize){
        if(mBuf[mBufDataLen - 1] != '\n'){
            mNoLineBreakAtEnd = true;
        }
    }
}

void FqReader::init(){
    if(util::endsWith(mFileName, ".gz")){
        mGzipFile = ::gzopen(mFileName.c_str(), "r");
        mZipped = true;
        ::gzrewind(mGzipFile);
    }else{
        if(mFileName == "/dev/stdin"){
            mFile = stdin;
        }else{
            mFile = std::fopen(mFileName.c_str(), "rb");
        }
        if(mFile == NULL){
            util::errorExit("Failed to open file: " + mFileName);
        }
        mZipped = false;
    }
    readToBuf();
}

void FqReader::getBytes(size_t& bytesRead, size_t& bytesTotal){
    if(mZipped){
        bytesRead = ::gzoffset(mGzipFile);
    }else{
        bytesRead = std::ftell(mFile);
    }
    
    // use another ifstream without affecting the current reader
    std::ifstream is(mFileName);
    is.seekg(0, is.end);
    bytesTotal = is.tellg();
}

void FqReader::clearLineBreaks(char* line){
    // trim \n, \r or \r\n in the tail
    int len = std::strlen(line);
    if(len >= 2){
        if(line[len - 1] == '\n' || line[len - 1] == '\r'){
            line[len - 1] = '\0';
            if(line[len - 2] == '\r'){
                line[len - 2] = '\0';
            }
        }
    }
}

std::string FqReader::getLine(){
    int start = mBufUsedLen;
    int end = start;
    
    // look for '\r' or '\n' until the end of mBuf
    while(end < mBufDataLen){
        if(mBuf[end] != '\r' && mBuf[end] != '\n'){
            ++end;
        }else{
            break;
        }
    }

    // if '\r' or '\n' found(this line well contained in this mBuf)
    // or this is the last mBuf of file
    if(end < mBufDataLen || mBufDataLen < mFqBufSize){
        int len = end - start;
        std::string line(mBuf + start, len);
        ++end;
        if(end < mBufDataLen - 1 && mBuf[end] == '\n'){
            ++end;
        }
        mBufUsedLen = end;
        return line;
    }

    // if '\r' or '\n' not found && this is not the last mBuf of file
    // then this line is not contained in this mBuf, we should read new mBuf
    std::string str(mBuf + start, mBufDataLen - start);
    while(true){
        readToBuf();
        start = 0;
        end = 0;
        // look for '\r' or '\n' until the end of mBuf
        while(end < mBufDataLen){
            if(mBuf[end] != '\r' && mBuf[end] != '\n'){
                ++end;
            }else{
                break;
            }
        }

        // if '\r' or '\n' found(this line well contained in this mBuf)
        // or this is the last mBuf of file
        if(end < mBufDataLen || mBufDataLen < mFqBufSize){
            int len = end - start;
            str.append(mBuf + start, len);
            ++end;
            if(end < mBufDataLen - 1 && mBuf[end] == '\n'){
                ++end;
            }
            mBufUsedLen = end;
            return str;
        }

        // if '\r' or '\n' not found && this is not the last mBuf of file
        // then this line is not contained in this mBuf, we should read new mBuf
        str.append(mBuf + start, mBufDataLen - start);
    }
    return std::string();
}

bool FqReader::eof(){
    if(mZipped){
        return ::gzeof(mGzipFile);
    }else{
        return std::feof(mFile);
    }
}

Read* FqReader::read(){
    if(mZipped && mGzipFile == NULL){
        return NULL;
    }
    if(mBufUsedLen >= mBufDataLen && eof()){
        return NULL;
    }

    std::string name = getLine();
    while((name.empty() && !(mBufUsedLen >= mBufDataLen && eof())) || (!name.empty() && name[0] !='@')){
        name = getLine();
    }
    if(name.empty()){
        return NULL;
    }

    std::string sequence = getLine();
    std::string strand = getLine();
    // some fq has no quality, then construct the quality string with all 'K'
    if(!mHasQuality){
        std::string quality = std::string(sequence.length(), 'K');
        return new Read(name, sequence, strand, quality, mPhread64);
    }else{
        std::string quality = getLine();
        if(quality.length() != sequence.length()){
            std::cerr << "Error: base sequnce and quality sequence have different length: \n";
            std::cerr << name << "\n";
            std::cerr << sequence << "\n";
            std::cerr << quality << "\n";
            std::cerr << strand << "\n";
            return NULL;
        }
        return new Read(name, sequence, strand, quality, mPhread64);
    }
    return NULL;
}

void FqReader::close(){
    if(mZipped && mGzipFile){
        ::gzclose(mGzipFile);
        mGzipFile = NULL;
    }else if(mFile){
        std::fclose(mFile);
        mFile = NULL;
    }else{
        return;
    }
}


bool FqReader::isZippedFq(const std::string& filename){
    // just use the suffix of filename to determine its format
    if(util::endsWith(filename, ".fastq.gz") || util::endsWith(filename, ".fq.gz") ||
       util::endsWith(filename, ".fasta.gz") || util::endsWith(filename, ".fa.gz")){
        return true;
    }
    return false;
}

bool FqReader::isfq(const std::string& filename){
    // just use the suffix of filename to determine its format
    if(util::endsWith(filename, ".fastq") || util::endsWith(filename, ".fq") ||
       util::endsWith(filename, ".fasta") || util::endsWith(filename, ".fa")){
        return true;
    }
    return false;
}

bool FqReader::isZipped(){
    return mZipped;
}

FqReaderPair::FqReaderPair(const std::string& lname, const std::string& rname, const bool& hasQual,
                           const bool& mPhread64, const bool& interleaved){
    mInterleaved = interleaved;
    left = new FqReader(lname, hasQual, mPhread64);
    if(interleaved){
        right = NULL;
    }else{
        right = new FqReader(rname, hasQual, mPhread64);
    }
}

FqReaderPair::~FqReaderPair(){
    if(left){
        delete left;
        left = NULL;
    }
    if(right){
        delete right;
        right = NULL;
    }
}

ReadPair* FqReaderPair::read(){
    Read* pr1 = left->read();
    Read* pr2 = NULL;
    if(mInterleaved){
        pr2 = left->read();
    }else{
        pr2 = right->read();
    }
    if(!pr1 || !pr2){
        return NULL;
    }else{
        return new ReadPair(pr1, pr2);
    }
}
