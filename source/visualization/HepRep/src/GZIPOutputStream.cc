// Copyright FreeHEP, 2005.

#include "cheprep/GZIPOutputStreamBuffer.h"
#include "cheprep/GZIPOutputStream.h"

/**
 * @author Mark Donszelmann
 * @version $Id: GZIPOutputStream.cc,v 1.3 2005-05-25 23:22:25 duns Exp $
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
