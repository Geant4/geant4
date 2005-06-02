// Copyright FreeHEP, 2005.
#ifndef ABSTRACTXMLWRITER_H
#define ABSTRACTXMLWRITER_H 1

#include "cheprep/config.h"

#include <string>

/**
 * @author Mark Donszelmann
 * @version $Id: AbstractXMLWriter.h,v 1.3 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {
    
    class AbstractXMLWriter {
    
        public:
            AbstractXMLWriter(std::string defaultNameSpace) : defaultNameSpace(defaultNameSpace) {
            }
            virtual ~AbstractXMLWriter() {
            }

            virtual void openTag(std::string ns, std::string name) = 0;
            virtual void printTag(std::string ns, std::string name) = 0;
            virtual void setAttribute(std::string ns, std::string name, std::string value) = 0;
            virtual void setAttribute(std::string ns, std::string name, double value) = 0;
                        
            virtual void close() = 0;
            virtual void openDoc(std::string version = "1.0", std::string encoding = "", bool standalone = false) = 0;
            virtual void closeDoc(bool force = false) = 0;
            virtual void openTag(std::string name) = 0;
            virtual void closeTag() = 0;
            virtual void printTag(std::string name) = 0;
            virtual void setAttribute(std::string name, char* value) = 0;
            virtual void setAttribute(std::string name, std::string value) = 0;
            virtual void setAttribute(std::string name, std::vector<double> value) = 0;
            virtual void setAttribute(std::string name, int64 value) = 0;
            virtual void setAttribute(std::string name, int value) = 0;
            virtual void setAttribute(std::string name, bool value) = 0;
            virtual void setAttribute(std::string name, double value) = 0;
        
        protected:
            std::string defaultNameSpace;
    };

} // cheprep

#endif  // ABSTRACTXMLWRITER_H
