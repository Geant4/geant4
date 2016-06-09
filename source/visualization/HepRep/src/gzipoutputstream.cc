// Copyright FreeHEP, 2005.

#include "cheprep/GZIPOutputStreamBuffer.h"
#include "cheprep/GZIPOutputStream.h"

/**
 * @author Mark Donszelmann
 * @version $Id: GZIPOutputStream.cc,v 1.4 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

    using namespace std;

    GZIPOutputStream::GZIPOutputStream(ostream &os)
                : std::ostream(NULL) {
      
        buffer = new GZIPOutputStreamBuffer(os.rdbuf()); 
        init(buffer);   
    }


    void GZIPOutputStream::setFilename(const string &filename) {
        buffer->setFilename(filename);
    }

    void GZIPOutputStream::setComment(const string &comment) {
        buffer->setComment(comment);
    }

    void GZIPOutputStream::close() {
        buffer->close();
    }


    GZIPOutputStream::~GZIPOutputStream() {
        delete buffer;
    }

} // cheprep
