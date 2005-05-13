#ifndef CHEPREP_GZIPOUTPUTSTREAMBUF_H
#define CHEPREP_GZIPOUTPUTSTREAMBUF_H

#include <string>

#include "cheprep/DeflateOutputStreamBuffer.h"

namespace cheprep {

    class GZIPOutputStreamBuffer : public DeflateOutputStreamBuffer {

        public:

            GZIPOutputStreamBuffer( std::streambuf *outbuf );

            void setFilename( const std::string &filename );
            void setComment( const std::string &comment );

            void close() ;

            virtual ~GZIPOutputStreamBuffer() ;

        protected:
            virtual int overflow( int c = EOF ) ;

        private:
            void writeHeader();
            void writeTrailer();
  
            std::string filename;
            std::string comment;
            bool open;
    };


} // cheprep

#endif // CHEPREP_GZIPOUTPUTSTREAMBUF_H
