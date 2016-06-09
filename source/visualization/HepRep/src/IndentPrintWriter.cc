// Copyright FreeHEP, 2005.

#include <fstream>

#include "cheprep/IndentPrintWriter.h"

using namespace std;

/**
 * @author Mark Donszelmann
 * @version $Id: IndentPrintWriter.cc,v 1.14 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

IndentPrintWriter::IndentPrintWriter(ostream* out, int level)
    : out(out), 
      closed(false), 
      indentLevel(level),
      indented(false),
      indentString("  ") {
}

IndentPrintWriter::~IndentPrintWriter() {
}

void IndentPrintWriter::close() {
    if (!closed) {
        out->flush();
        closed = true;
    }
}

IndentPrintWriter& IndentPrintWriter::operator<< (const std::string & s) {
    if (!indented) doIndent();
    *out << s;
    return *this;
}

IndentPrintWriter& IndentPrintWriter::operator<< (ostream& (*)(ostream&)) {
    *out << endl;
    indented = false;
    return *this;
}

void IndentPrintWriter::println(const string & s) {
	*this << s << endl;
}

void IndentPrintWriter::print(const string & s) {
    *this << s;
}

void IndentPrintWriter::println() {
    *out << endl;
	indented = false;
}

void IndentPrintWriter::doIndent() {
	for (int i=0; i<indentLevel; i++) {
	    *out << indentString;
	}
	indented = true;
}

void IndentPrintWriter::indent() {
    indentLevel++;
}

void IndentPrintWriter::outdent() {
    indentLevel--;
}

int IndentPrintWriter::getIndent() const {
    return indentLevel;
}

void IndentPrintWriter::setIndent(const int level) {
    indentLevel = level;
}

string IndentPrintWriter::getIndentString() const {
    return indentString;
}

void IndentPrintWriter::setIndentString(const string & indent) {
    indentString = indent;
}

} // cheprep
