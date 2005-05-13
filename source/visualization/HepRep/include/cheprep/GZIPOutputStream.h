#ifndef CHEPREP_GZIPOUTPUTSTREAM_H
#define CHEPREP_GZIPOUTPUTSTREAM_H

#include <string>

#include "cheprep/GZIPOutputStreamBuffer.h"

namespace cheprep {

    class GZIPOutputStream : public std::ostream {
        
        public:

            GZIPOutputStream(std::ostream &os);

            void setFilename(const std::string &filename);
            void setComment(const std::string &comment);
  
            void close();

            virtual ~GZIPOutputStream();

        private:
            GZIPOutputStreamBuffer *buffer;
    };
 
} // cheprep.

#endif // CHEPREP_GZIPOUTPUTSTREAM_H
