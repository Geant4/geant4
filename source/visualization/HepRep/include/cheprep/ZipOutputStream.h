// Copyright FreeHEP, 2005.
#ifndef CHEPREP_ZIPOUTPUTSTREAM_H
#define CHEPREP_ZIPOUTPUTSTREAM_H

#include <string>
#include <iostream>
#include <vector>


/**
 * @author Mark Donszelmann
 * @version $Id: ZipOutputStream.h 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

    class ZipOutputStreamBuffer;

    class ZipOutputStream : public std::ostream {

        public:

            ZipOutputStream(std::ostream& os);
  
            void closeEntry();

            void close();

            void putNextEntry(const std::string& name, bool compress);
            
            void setComment(const std::string& comment);

            virtual ~ZipOutputStream();

        private:
            ZipOutputStreamBuffer* buffer;
    };
 
} // cheprep

#endif // CHEPREP_ZIPOUTPUTSTREAM_H
