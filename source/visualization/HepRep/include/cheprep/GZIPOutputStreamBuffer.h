// Copyright FreeHEP, 2005.
#ifndef CHEPREP_GZIPOUTPUTSTREAMBUF_H
#define CHEPREP_GZIPOUTPUTSTREAMBUF_H

#include <string>

#include "cheprep/DeflateOutputStreamBuffer.h"

/**
 * @author Mark Donszelmann
 * @version $Id: GZIPOutputStreamBuffer.h 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

    class GZIPOutputStreamBuffer : public DeflateOutputStreamBuffer {

        public:

            GZIPOutputStreamBuffer( std::streambuf *outbuf );

            int overflow(int);

            void setFilename( const std::string &filename );
            void setComment( const std::string &comment );

            void close() ;

            virtual ~GZIPOutputStreamBuffer() ;

        private:
            void writeHeader();
            void writeTrailer();
  
            std::string filename;
            std::string comment;
            bool open;
    };


} // cheprep

#endif // CHEPREP_GZIPOUTPUTSTREAMBUF_H
