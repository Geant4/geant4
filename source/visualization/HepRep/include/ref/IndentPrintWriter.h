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
 */
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

#endif

