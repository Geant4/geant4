// Copyright FreeHEP, 2005.
#ifndef CHEPREP_ZIPOUTPUTSTREAMBUFFER_H
#define CHEPREP_ZIPOUTPUTSTREAMBUFFER_H

#include <string>
#include <iostream>
#include <vector>

#include "cheprep/DeflateOutputStreamBuffer.h"

/**
 * @author Mark Donszelmann
 * @version $Id: ZipOutputStreamBuffer.h,v 1.9 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

    class ZipEntry;

    class ZipOutputStreamBuffer : public DeflateOutputStreamBuffer {

        public:

            ZipOutputStreamBuffer(std::streambuf* buffer);
  
            int overflow(int c);
    
            void closeEntry();

            void close();

            void putNextEntry(const std::string& name, bool compress);

            void setMethod(int method);
            
            void setComment(const std::string& comment);

            virtual ~ZipOutputStreamBuffer();

        private:
            std::string comment;

            bool closed;
            ZipEntry* entry;
            std::vector<ZipEntry*>* entries;
                                
            static const unsigned int LOCSIG = 0x04034b50;
            static const unsigned int EXTSIG = 0x08074b50;
            static const unsigned int CENSIG = 0x02014b50;
            static const unsigned int ENDSIG = 0x06054b50;
        
            static const unsigned int VERSIONMADE           = 0x0014;
            static const unsigned int VERSIONEXTRACT        = 0x0014;
            static const unsigned int GENFLAG               = 0x0008;
    };
 
} // cheprep

#endif // CHEPREP_ZIPOUTPUTSTREAMBUFFER_H
