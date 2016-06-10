// Copyright FreeHEP, 2005.
#ifndef CHEPREP_INDENTPRINTWRITER_H
#define CHEPREP_INDENTPRINTWRITER_H 1

#include "cheprep/config.h"

#include <iostream>
#include <string>

/**
 * A PrintWriter that keeps track of an indentation level
 * and indents the output appropriately.
 *
 * <b>Warning:</b> Only print and println methods taking strings have been overriden,
 * print, println methods taking other arguments may not be indented properly.
 *
 * @author Mark Donszelmann
 * @version $Id: IndentPrintWriter.h 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

class IndentPrintWriter {

	public:
	    IndentPrintWriter(std::ostream* out, int level = 0);
        virtual ~IndentPrintWriter();

        void close();
        IndentPrintWriter& operator<< (const std::string & s);
        IndentPrintWriter& operator<< (std::ostream& (*pf)(std::ostream&));
	    void println(const std::string & s);
        void print(const std::string & s);
	    void println();
	    void indent();
	    void outdent();
	    int getIndent() const;
        void setIndent(const int level);
        std::string getIndentString() const;
        void setIndentString(const std::string & indentString);

    private:
        void doIndent();

        std::ostream* out;
        bool closed;
        int indentLevel;
        bool indented;
        std::string indentString;
};

} // cheprep


#endif

