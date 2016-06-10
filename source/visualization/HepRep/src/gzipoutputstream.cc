// Copyright FreeHEP, 2005.

#include "cheprep/GZIPOutputStreamBuffer.h"
#include "cheprep/GZIPOutputStream.h"

/**
 * @author Mark Donszelmann
 * @version $Id: GZIPOutputStream.cc 66373 2012-12-18 09:41:34Z gcosmo $
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
