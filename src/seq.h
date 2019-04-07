#ifndef SEQ_H
#define SEQ_H

#include <cstdio>
#include <string>
#include <cstdlib>
#include <iostream>

/** class to represent a nucleotide sequence */
class Seq{
    public:
        /** default constructor */
        Seq() = default;
        
        /** construct a Seq object with a sequence s */
        Seq(const std::string& s) : seqStr(s){ }
        
        /** get the length of this Seq object */
        inline int length(){
            return seqStr.length();
        }
        
        /** get the reverse complementary Seq object */
        inline Seq reverseComplement(){
            std::string rstr(seqStr.length(), 0);
            char chc = '\0';
            for(size_t i = 0; i < seqStr.length(); ++i){
                switch(seqStr[i]){
                    case 'A': case 'a':
                        chc = 'T';
                        break;
                    case 'T': case 't':
                        chc = 'A';
                        break;
                    case 'C': case 'c':
                        chc = 'G';
                        break;
                    case 'G': case 'g':
                        chc = 'C';
                        break;
                    default:
                        chc = 'N';
                        break;
                }
                rstr[rstr.length() - 1 - i] = chc;
            }
            return Seq(rstr);
        }
        
        /** output Seq to ostream
         * @param os std::ostream object
         * @param s Seq object
         */
        friend inline std::ostream& operator<<(std::ostream& os, const Seq& s){
            os << s.seqStr;
            return os;
        }
        
        /** operator to get the reverse complement Seq object */
        inline Seq operator~(){
            return this->reverseComplement();
        }

    public:
        std::string seqStr; ///< store the nucleotide sequence
};

#endif
