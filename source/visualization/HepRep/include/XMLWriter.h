#ifndef XMLWRITER_H
#define XMLWRITER_H 1

#include "FreeHepTypes.h"

#include <iostream>
#include <map>
#include <stack>
#include <vector>
#include <string>

#include "IndentPrintWriter.h"

/**
 * A class that makes it easy to write XML documents.
 *
 * @author Tony Johnson
 * @author Mark Donszelmann
 * @version $Id: XMLWriter.h,v 1.7 2002-11-19 21:54:35 duns Exp $
 */
class XMLWriter {

    public:
        XMLWriter(std::ostream* out, std::string indentString = "  ", std::string defaultNameSpace = "");
        virtual ~XMLWriter();
        void close();
        void openDoc(std::string version = "1.0", std::string encoding = "", bool standalone = false);
        void referToDTD(std::string name, std::string pid, std::string ref);
        void referToDTD(std::string name, std::string system);
        void closeDoc();
        void printComment(std::string comment);
        void print(std::string text);
        void println(std::string text);
        void openTag(std::string ns, std::string name);
        void openTag(std::string name);
        void closeTag();
        void printTag(std::string ns, std::string name);
        void printTag(std::string name);
        void setAttribute(std::string name, std::string value);
        void setAttribute(std::string ns, std::string name, std::string value);
        void setAttribute(std::string name, double value);
        void setAttribute(std::string ns, std::string name, double value);
        void printAttributes(int tagLength);
        std::string normalize(std::string s);
        std::string normalizeText(std::string s);
        void checkNameValid(std::string s);

    protected:
        bool closed;
	    IndentPrintWriter* writer;
        std::string defaultNameSpace;

    private:
        std::string* dtdName;
	    std::map<std::string, std::string> attributes;
	    std::stack<std::string> openTags;
};

#endif
