
#include "cheprep/GZIPOutputStreamBuffer.h"
#include "cheprep/GZIPOutputStream.h"

namespace cheprep {

    using namespace std;

    GZIPOutputStream::GZIPOutputStream(ostream &os)
                : std::ostream(NULL) {
      
        buffer = new GZIPOutputStreamBuffer(os.rdbuf());    
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
