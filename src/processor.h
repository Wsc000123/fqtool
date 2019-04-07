#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <string>
#include "options.h"
#include "seprocessor.h"
#include "peprocessor.h"

/** class to do se/pe processing */
class Processor{
    public:
        /* construct a Processor object 
         * @param opt pointer to Options object
         */
        Processor(Options* opt);
        
        /* destroy a Processor object */
        ~Processor();

        /** process fastq */
        bool process();
    public:
        Options* mOptions;  ///< pointer to Options object
};

#endif
