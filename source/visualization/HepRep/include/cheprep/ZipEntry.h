// Copyright FreeHEP, 2005.
#ifndef CHEPREP_ZIPENTRY_H
#define CHEPREP_ZIPENTRY_H

#include <string>
#include <vector>

/**
 * @author Mark Donszelmann
 * @version $Id: ZipEntry.h 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

    class ZipEntry {

        public:

            std::string name;
            int method;
            int date;
            int time;
            unsigned int crc;
            unsigned int size;
            unsigned int csize;
            long offset;
            long data;
            
            ZipEntry() {
            }
  
            virtual ~ZipEntry() {
            }
    };
 
} // cheprep

#endif // CHEPREP_ZIPENTRY_H
