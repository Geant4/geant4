
#include <fstream>

#include "IndentPrintWriter.h"

using namespace std;

IndentPrintWriter::IndentPrintWriter(ostream* out, int level)
    : out(out), indentLevel(level) {

    indentString = "  ";
    indented = false;
}

IndentPrintWriter::~IndentPrintWriter() {
}

void IndentPrintWriter::close() {
    out->flush();
    ofstream* fout = dynamic_cast<ofstream *>(out);
    if (fout != NULL) {
        fout->close();
    }
}

IndentPrintWriter& IndentPrintWriter::operator<< (const char* s) {
    if (!indented) doIndent();
    *out << s;
    return *this;
}

IndentPrintWriter& IndentPrintWriter::operator<< (ostream& (*pf)(ostream&)) {
    *out << endl;
    indented = false;
    return *this;
}

void IndentPrintWriter::println(string s) {
	*this << s.c_str() << endl;
}

void IndentPrintWriter::print(string s) {
    *this << s.c_str();
}

void IndentPrintWriter::println() {
    *out << endl;
	indented = false;
}

void IndentPrintWriter::doIndent() {
	for (int i=0; i<indentLevel; i++) {
	    *out << indentString.c_str();
	}
	indented = true;
}

void IndentPrintWriter::indent() {
    indentLevel++;
}

void IndentPrintWriter::outdent() {
    indentLevel--;
}

int IndentPrintWriter::getIndent() {
    return indentLevel;
}

void IndentPrintWriter::setIndent(int level) {
    indentLevel = level;
}

string IndentPrintWriter::getIndentString() {
    return indentString;
}

void IndentPrintWriter::setIndentString(string indent) {
    indentString = indent;
}

