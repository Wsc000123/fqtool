#ifndef UMI_PROCESSOR_H
#define UMI_PROCESSOR_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include "options.h"
#include "read.h"

/** Class to process umi */
class UmiProcessor{
    public:
        /** construct a UmiProcessor */
        UmiProcessor(Options* opt);

        /** destruct a UmiProcessor */
        ~UmiProcessor();

        /** find the umi according position arguments provided\n
         * and add umi (and qual) raw tag to read name
         * @param r1 pointer to Read (read1)
         * @param r2 pointer to Read (read2)
         */
        void process(Read* r1, Read* r2 = NULL);
        
        /** add umi string to read 
         * @param r pointer to Read
         * @param tag tag string
         */ 
        void addTagToName(Read* r, const std::string& tag);
    private:
        Options* mOptions;
        static constexpr int UMI_LOC_INDEX1 = 1;    ///< UMI is first index of read1
        static constexpr int UMI_LOC_INDEX2 = 2;    ///< UMI is first index of read2
        static constexpr int UMI_LOC_READ1 = 3;     ///< UMI is the first few bases of read1
        static constexpr int UMI_LOC_READ2 = 4;     ///< UMI is the first few bases of read2
        static constexpr int UMI_LOC_PER_INDEX = 5; ///< UMI is the first index of each read
        static constexpr int UMI_LOC_PER_READ = 6;  ///< UMI is the first few bases of each read
};

#endif
