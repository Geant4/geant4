#ifndef INDENTPRINTWRITER_H
#define INDENTPRINTWRITER_H 1

#include "FreeHepTypes.h"

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
 * @version $Id: IndentPrintWriter.h,v 1.7 2002-11-19 21:53:57 duns Exp $
 */
class IndentPrintWriter {

	public:
	    IndentPrintWriter(std::ostream* out, int level = 0);
        virtual ~IndentPrintWriter();

        void close();
        IndentPrintWriter& operator<< (const char *s);
        IndentPrintWriter& operator<< (std::ostream& (*pf)(std::ostream&));
	    void println(std::string s);
        void print(std::string s);
	    void println();
	    void indent();
	    void outdent();
	    int getIndent();
        void setIndent(int level);
        std::string getIndentString();
        void setIndentString(std::string indentString);

    private:
        void doIndent();

	std::ostream* out;
        int indentLevel;
	bool indented;
	std::string indentString;
};

#endif

