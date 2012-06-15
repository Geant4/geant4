// Copyright FreeHEP, 2005.

#include "cheprep/GZIPOutputStreamBuffer.h"

/**
 * @author Mark Donszelmann
 * @version $Id: GZIPOutputStreamBuffer.cc,v 1.4 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

    using namespace std;

    GZIPOutputStreamBuffer::GZIPOutputStreamBuffer( streambuf *aBuffer)
        : DeflateOutputStreamBuffer(aBuffer),
          open(false) {
            
        init(true);
    }
    
    void GZIPOutputStreamBuffer::setFilename( const string &name ) {
        filename = name ;
    }
    
    void GZIPOutputStreamBuffer::setComment( const string &c ) {
        comment = c ;
    }
    
    void GZIPOutputStreamBuffer::close() {
        if (!open) return;

        finish();         
        writeTrailer();
      
        open = false ;
    }
    
    GZIPOutputStreamBuffer::~GZIPOutputStreamBuffer() {
        close() ;
    }
    
    int GZIPOutputStreamBuffer::overflow( int c ) {
        if (!open) {
            writeHeader();
            open = true;
        }
        return DeflateOutputStreamBuffer::overflow( c ) ;
    }
        
    void GZIPOutputStreamBuffer::writeHeader() {
        unsigned char flag = 0x00;
        flag |= (filename == "") ? 0x00 : 0x08;
        flag |= (comment  == "") ? 0x00 : 0x10;
    
        putUB(0x1f);  // Magic #
        putUB(0x8b);  // Magic #        
        putUB(0x08);  // Deflater.DEFLATED
        putUB(flag);  // FLG
        putUI(0x00000000);  // MTIME
        putUB(0x00);  // XFLG
        putUB(0x00);  // OS
            
        if (filename != "") {
            putS(filename); // Filename
            putUB(0x00);
        }
        
        if (comment != "") {
            putS(comment); // Comment
            putUB(0x00);
        }
    }
    
    void GZIPOutputStreamBuffer::writeTrailer() {
        putUI(getCRC());
        putUI(getSize());
    }
    
} // cheprep
