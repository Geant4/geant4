// Copyright FreeHEP, 2005.

#include <iostream>
#include <ctime>
#include <vector>

#include "cheprep/ZipOutputStream.h"
#include "cheprep/ZipOutputStreamBuffer.h"
#include "cheprep/ZipEntry.h"

/**
 * @author Mark Donszelmann
 * @version $Id: ZipOutputStreamBuffer.cc 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

    ZipOutputStreamBuffer::ZipOutputStreamBuffer(std::streambuf* aBuffer)
        : DeflateOutputStreamBuffer(aBuffer), 
          comment("") {
        
        closed = false;
        entry = NULL;
        
        entries = new std::vector<ZipEntry*>();
    }

    int ZipOutputStreamBuffer::overflow(int c) {
        return DeflateOutputStreamBuffer::overflow(c);
    }
            
    void ZipOutputStreamBuffer::closeEntry() {
        if (closed) return;
        if (entry == NULL) return;
        
        finish();
                
        entry->crc = getCRC();
        // NOTE: pos()::operator- is ambiguous on gcc 4.0, so convert to long first
        entry->csize = (long)pos() - entry->data;
        entry->size = getSize();
        putUI(EXTSIG);
        putUI(entry->crc);            // crc
        putUI(entry->csize);          // compressed size
        putUI(entry->size);           // uncompressed size
        entry = NULL;
    }


    void ZipOutputStreamBuffer::close() {
        if (closed) return;
        closeEntry();
        
        long dirStart = pos();
        for (std::vector<ZipEntry*>::iterator i=entries->begin(); i != entries->end(); i++) {
            entry = *i;
            putUI(CENSIG);
            putUS(VERSIONMADE);                   // made by version
            putUS(VERSIONEXTRACT);                // extraction version
            putUS(GENFLAG);                       // general purpose bit flag
            putUS(entry->method);                 // compression method
            putUS(entry->time);                   // time
            putUS(entry->date);                   // date
            putUI(entry->crc);                    // crc
            putUI(entry->csize);                  // compressed size
            putUI(entry->size);                   // uncompressed size
            putUS(entry->name.length());          // file name length
            putUS(0);                             // extra field length
            putUS(0);                             // file comment length
            putUS(0);                             // disk number start
            putUS(0);                             // internal file attributes
            putUI(0);                             // external file attributes
            putUI(entry->offset);                 // relative offset of local file header
            putS(entry->name);                    // file name
            
            // delete entry
            delete entry;
            entry = NULL;
        }
        // NOTE: pos()::operator- is ambiguous on gcc 4.0, so convert to long first
        long dirSize = (long)pos() - dirStart;
        putUI(ENDSIG);
        putUS(0);                             // number of this disk
        putUS(0);                             // number of the disk with the start of central directory
        putUS(entries->size());               // total number of entries in the central directory of this disk
        putUS(entries->size());               // total number of entries in the central directory
        putUI(dirSize);                       // size of the central directory
        putUI(dirStart);                      // offset of the start of central directory with respect to the starting disk number
        putUS(comment.length());              // .ZIP file comment length
        putS(comment);                        // .ZIP file comment
        
        delete entries;
        entries = NULL;
        closed = true;
    }

    void ZipOutputStreamBuffer::putNextEntry(const std::string& name, bool compress) {
        if (closed) return;
        
        closeEntry();

#ifdef CHEPREP_NO_ZLIB
        compress = false;
#endif
        init(compress);

        entry = new ZipEntry();
        entries->push_back(entry);
        
        entry->name = name;
        entry->method = compress ? 8 : 0;
        
        time_t ltime;
        time( &ltime );
        struct tm *utc = gmtime( &ltime );
        entry->date = (utc->tm_year - 80) << 9 | (utc->tm_mon + 1) << 5 | utc->tm_mday;
        entry->time = utc->tm_hour << 11 | utc->tm_min << 5 | utc->tm_sec >> 1;

        entry->offset = (long)pos();
        putUI(LOCSIG);                        // signature
        putUS(VERSIONEXTRACT);                // extraction version
        putUS(GENFLAG);                       // general purpose bit flag
        putUS(entry->method);                 // compression method
        putUS(entry->time);                   // time
        putUS(entry->date);                   // date
        putUI(0x00000000);                    // crc ATEND
        putUI(0x00000000);                    // compressed size ATEND
        putUI(0x00000000);                    // uncompressed size ATEND
        putUS(entry->name.length());          // file name length
        putUS(0);                             // extra field length
        putS(entry->name);                    // file name
        entry->data = (long)pos();
        entry->crc = 0;
    }
            
    void ZipOutputStreamBuffer::setComment(const std::string& c) {
        if (closed) return;
        comment = c;        
    }

    ZipOutputStreamBuffer::~ZipOutputStreamBuffer() {
        close();
    }

} // cheprep
