// Copyright FreeHEP, 2005.
#ifndef CHEPREP_XMLWRITER_H
#define CHEPREP_XMLWRITER_H 1

#include "cheprep/config.h"

#include <iostream>
#include <map>
#include <stack>
#include <vector>
#include <string>

#include "cheprep/AbstractXMLWriter.h"
#include "cheprep/IndentPrintWriter.h"

/**
 * @author Mark Donszelmann
 * @version $Id: XMLWriter.h 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

class XMLWriter : public AbstractXMLWriter {

    public:
        XMLWriter(std::ostream* out, std::string indentString = "  ", std::string defaultNameSpace = "");
        virtual ~XMLWriter();
        void close();
        void openDoc(std::string version = "1.0", std::string encoding = "", bool standalone = false);
        void referToDTD(std::string name, std::string pid, std::string ref);
        void referToDTD(std::string name, std::string system);
        void closeDoc(bool force = false);
        void printComment(std::string comment);
        void printPlain(std::string text);
        void print(std::string text);
        void println(std::string text);
        void openTag(std::string name);
        void closeTag();
        void printTag(std::string name);
        void setAttribute(std::string name, char* value);
        void setAttribute(std::string name, std::string value);
        void setAttribute(std::string name, std::vector<double> value);
        void setAttribute(std::string name, int64 value);
        void setAttribute(std::string name, int value);
        void setAttribute(std::string name, bool value);
        void setAttribute(std::string name, double value);
        void printAttributes(int tagLength);
        std::string normalize(std::string s);
        std::string normalizeText(std::string s);
        void checkNameValid(std::string s);

        //
        // Can be removed when we can properly inherit those (since names are equal to overloaded ones).
        // 
        void openTag(std::string ns, std::string name) {
            openTag(ns == defaultNameSpace ? name : ns.append(":").append(name));
        }
        void printTag(std::string ns, std::string name) {
            printTag(ns == defaultNameSpace ? name : ns.append(":").append(name));
        }
        void setAttribute(std::string ns, std::string name, std::string value) {
            setAttribute(ns.append(":").append(name), value);
        }
        void setAttribute(std::string ns, std::string name, double value) {
            setAttribute(ns.append(":").append(name), value);
        }



    protected:
        bool closed;
	    IndentPrintWriter* writer;

    private:
        std::string dtdName;
	    std::map<std::string, std::string> attributes;
	    std::stack<std::string> openTags;
};

} // cheprep

#endif
