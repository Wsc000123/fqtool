#ifndef READ_H
#define READ_H

#include "seq.h"
#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

/** Class to represent a ngs read record */
class Read{
    public:
        std::string name;   ///< read name
        Seq seq;            ///< read nucleotide sequence
        std::string strand; ///< read strand
        std::string quality;///< read quality sequence
        bool hasQuality;    ///< read has quality sequence if true

    public:
        /** default constructor of Read
         */
        Read() = default;

        /** Read constructor
         * @param rname read name
         * @param rseq read sequence
         * @param rstrand read strand
         * @param rqual read quality sequence
         * @param phread64 quality is encoded in phread64 or not
         */
        Read(const std::string& rname, const std::string& rseq, const std::string& rstrand, const std::string& rqual, const bool& phread64 = false) :
            name(rname), seq(rseq), strand(rstrand), quality(rqual), hasQuality(true){
                if(phread64){
                    convertPhread64To33();
                }
            }
        
        /** Read constructor
         * @param rname read name
         * @param oseq Seq object
         * @param rstrand read strand
         * @param rqual read quality sequence
         * @param phread64 quality is encoded in phread64 or not
         */
        Read(const std::string& rname, const Seq& oseq, const std::string& rstrand, const std::string& rqual, const bool& phread64 = false) :
            name(rname), seq(oseq), strand(rstrand), quality(rqual), hasQuality(true){
                if(phread64){
                    convertPhread64To33();
                }
            }

        /** Read constructor
         * @param rname read name
         * @param oseq Seq object
         * @param rstrand read strand
         */
        Read(const std::string& rname, const Seq& oseq, const std::string& rstrand) :
            name(rname), seq(oseq), strand(rstrand), hasQuality(false){ }
        
        /** Read constructor
         * @param r Read object
         */
        Read(const Read& r) : name(r.name), seq(r.seq), strand(r.strand), quality(r.quality), hasQuality(r.hasQuality){}
        
        /** Read destructor */
        ~Read(){};

        /** convert quality in phread64 based encoding into phread33 based encoding */
        inline void convertPhread64To33(){
            for(size_t i = 0; i < quality.length(); ++i){
                quality[i] = std::max(33, quality[i] - (64 - 33));
            }
        }

        /** output a Read object to std::ostream
         * @param os an std::ostream
         * @param r a Read object
         * @return os
         */
        friend std::ostream& operator<<(std::ostream& os, const Read& r){
            os << r.name << "\n";
            os << r.seq.seqStr << "\n";
            os << r.strand << "\n";
            if(r.hasQuality){
                os << r.quality << "\n";
            }
            return os;
        }

        /** get a reverse complementary Read
         * @return pointer to the complementary Read
         */
        inline Read* reverseComplement(){
            Seq rseq = ~seq;
            std::string rqual(quality.rbegin(), quality.rend());
            std::string rstrand = (strand == "-" ? "+" : "-");
            return new Read(name, rseq, rstrand, rqual);
        }

        /** readname example '\@A00403:136:HFMYWDSXX:2:1101:7672:1000 1:N:0:GAGAGGCA+GAGAGGC'
         * get the first index of a Read in the Read name
         * @return first index of Read with two index or just the index of Read
         */
        std::string firstIndex(){
            int len = name.length();
            int end = len;
            // too short to hold an index
            if(len < 5){
                return "";
            }
            // index is at least 2 base long
            for(int i = len - 3; i >=0; --i){
                if(name[i] == '+'){
                    end = i - 1;
                }
                if(name[i] == ':'){
                    return name.substr(i + 1, end - i);
                }
            }
            return "";
        }

        /** get the last index of a Read in the Read name
         * @return last index of Read with two index or just the index of Read
         */
        inline std::string lastIndex(){
            int len = name.length();
            if(len < 5){
                return "";
            }
            for(int i = len - 3; i >=0; --i){
                //: for single end sequence; + for pair end sequence
                if(name[i] == ':' || name[i] == '+'){
                    return name.substr(i + 1, len - i);
                }
            }
            return "";
        }

        /** count bases in a Read with quality lower than a threshold (0-based)
         * @param lowQual threshold lower than it will be counted as low quality
         * @return number of bases with quality lower than lowQual
         */
        inline int lowQualCount(const int& lowQual=20){
            int count = 0;
            for(size_t i = 0; i < quality.length(); ++i){
                if(quality[i] < lowQual + 33){
                    ++count;
                }
            }
            return count;
        }

        /** get the length of a Read
         * @return the length of the Read
         */
        inline int length(){
            return seq.length();
        }

        /** convert a Read object to a string
         * @return a string representation of a Read
         */
        inline std::string toString(){
            return name + "\n" + seq.seqStr + "\n" + strand + "\n" + quality + "\n";
        }
        
        /** resize a Read to specified length
         * @param len length
         */
        inline void resize(int len){
            if(len > this->length() || len < 0){
                return;
            }
            seq.seqStr.resize(len);
            quality.resize(len);
        }

        /** trim a Read from front (5')
         * @param len length to be trimmed
         */
        inline void trimFront(int len){
            int orilen = length();
            len = std::min(len, orilen - 1);
            seq.seqStr = seq.seqStr.substr(len, orilen - len);
            quality = quality.substr(len, orilen - len);
        }
};

/** class to represent a pair of read */
class ReadPair{
    public:
        Read* left; ///< Pointer to read1
        Read* right;///< Pointer to read2

    public:
        /** Default ReadPair Constructor */
        ReadPair() = default;
        
        /** Construct a ReadPair from read1/2 pointers */
        ReadPair(Read* rLeft, Read* rRight) : left(rLeft), right(rRight){}
        
        /** Destroy a ReadPair object */
        ~ReadPair(){
            if(left){
                delete left;
                left = NULL;
            }
            if(right){
                delete right;
                right = NULL;
            }
        }

        /** merge a pair of Reads, ignore false indel caused by sequence error
         * @return merged Read
         */
        Read* merge();
};

#endif
