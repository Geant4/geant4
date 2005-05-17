// Copyright FreeHEP, 2005.
#ifndef CHEPREP_ZIPENTRY_H
#define CHEPREP_ZIPENTRY_H

#include <string>
#include <vector>

/**
 * @author Mark Donszelmann
 * @version $Id: ZipEntry.h,v 1.6 2005-05-17 22:07:07 duns Exp $
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
