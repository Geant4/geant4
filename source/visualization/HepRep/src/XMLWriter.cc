// Copyright FreeHEP, 2005.

#include "cheprep/config.h"

#include <cstdio>

#include "cheprep/DefaultHepRepAttValue.h"
#include "cheprep/XMLWriter.h"

using namespace std;

/**
 * @author Mark Donszelmann
 * @version $Id: XMLWriter.cc,v 1.12 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

XMLWriter::XMLWriter(ostream* out, string indentString, string aDefaultNameSpace)
    : AbstractXMLWriter(aDefaultNameSpace) {
    writer = new IndentPrintWriter(out);
    writer->setIndentString(indentString);
    closed = false;
    dtdName = "";
}

XMLWriter::~XMLWriter() {
    writer->close();
    delete writer;
}

void XMLWriter::close() {
    closeDoc();
    writer->close();
}

void XMLWriter::openDoc(string version, string encoding, bool standalone) {
    string indentString = writer->getIndentString();
    writer->setIndentString(indentString);

//    if (!XMLCharacterProperties.validVersionNum(version)) throw new RuntimeException("Invalid version number: "+version);
    *writer << "<?xml version=\"" << version.c_str() << "\" ";
    if (encoding.compare("") != 0) {
//        if (!XMLCharacterProperties.validEncName(encoding)) throw new RuntimeException("Invalid encoding name: "+encoding);
        *writer << "encoding=\"" << encoding.c_str() << "\" ";
    }
    if (standalone) {
        *writer << "standalone=\"yes\" ";
    }
    *writer << "?>";
    *writer << endl;
    writer->setIndentString(indentString);
}

void XMLWriter::referToDTD(string name, string pid, string ref) {
    if (dtdName != "") {
        cerr << "XMLWriter::ReferToDTD cannot be called twice" << endl;
    }
    dtdName = name;
    *writer << "<!DOCTYPE " << name.c_str() << " PUBLIC \"" << pid.c_str() << "\" \"" << ref.c_str() << "\">" << endl;
}

void XMLWriter::referToDTD(string name, string system) {
    if (dtdName != "") {
        cerr << "XMLWriter::ReferToDTD cannot be called twice";
    }
    dtdName = name;
    *writer << "<!DOCTYPE " << name.c_str() << " SYSTEM \"" << system.c_str() << "\">" << endl;
}

void XMLWriter::closeDoc(bool force) {
    if (!closed) {
        if (!openTags.empty()) {
            if (!force) cerr << "Not all tags were closed before closing XML document:" << endl;
            while (!openTags.empty()) {
                if (force) {
                    closeTag();
                } else {
                    cerr << "   </" << openTags.top().c_str() << ">" << endl;
                    openTags.pop();
                }
            }
        }
        closed = true;
    }
}

void XMLWriter::printComment(string comment) {
    if (comment.find("--") != string::npos) {
        cerr << "XMLWriter::printComment '--' sequence not allowed in comment" << endl;
    }
    *writer << "<!--" << normalizeText(comment).c_str() << "-->" << endl;
}

void XMLWriter::printPlain(string text) {
    *writer << text.c_str();
}

void XMLWriter::print(string text) {
    *writer << normalizeText(text).c_str();
}

void XMLWriter::println(string text) {
    print(text);
    *writer << endl;
}

void XMLWriter::openTag(string name) {
    checkNameValid(name);
    if (openTags.empty() && dtdName.compare("") && dtdName.compare(name)) {
        cerr << "XMLWriter::openTag(), First tag: '" << name << "' not equal to DTD id: '" << dtdName << "'" << endl;
    }
    *writer << "<" << name.c_str();
    printAttributes(name.length());
    *writer << ">" << endl;
    writer->indent();
    openTags.push(name);
}

void XMLWriter::closeTag() {
    if (openTags.empty()) {
        writer->close();
        cerr << "XMLWriter::closeTag(), No open tags" << endl;
    }
    string name = openTags.top();
    openTags.pop();
    writer->outdent();
    *writer << "</" << name.c_str() << ">" << endl;
}

void XMLWriter::printTag(string name) {
    checkNameValid(name);
    *writer << "<" << name.c_str();
    printAttributes(name.length());
    *writer << "/>" << endl;
}

void XMLWriter::setAttribute(string name, char* value) {
    setAttribute(name, (string)value);
}

void XMLWriter::setAttribute(string name, string value) {
    attributes[name] = value;
    // NOTE: never set type here
}

void XMLWriter::setAttribute(std::string name, std::vector<double> value) {
    if (name == "value") setAttribute("type", (std::string)"Color");
    setAttribute(name, DefaultHepRepAttValue::getAsString(value));
}

void XMLWriter::setAttribute(std::string name, int64 value) {
    if (name == "value") setAttribute("type", (std::string)"long");
    setAttribute(name, DefaultHepRepAttValue::getAsString(value));
}

void XMLWriter::setAttribute(std::string name, int value) {
    if (name == "showlabel") {
        string label = DefaultHepRepAttValue::toShowLabel(value);
        setAttribute("showlabel", label);
    } else {
        if (name == "value") setAttribute("type", (std::string)"int");
        setAttribute(name, DefaultHepRepAttValue::getAsString(value));
    }
}

void XMLWriter::setAttribute(std::string name, bool value) {
    if (name == "value") setAttribute("type", (std::string)"boolean");
    setAttribute(name, DefaultHepRepAttValue::getAsString(value));
}

void XMLWriter::setAttribute(string name, double value) {
    if (name == "value") setAttribute("type", (std::string)"double");
    setAttribute(name, DefaultHepRepAttValue::getAsString(value));
}

void XMLWriter::printAttributes(int tagLength) {
	int width = tagLength + 1;
	bool extraIndent = false;
	for (map<string,string>::iterator i = attributes.begin(); i != attributes.end(); i++) {
		string key = i->first;
		checkNameValid(key);
		string value = normalize(i->second);
		int length = key.length() + value.length() + 3;
		if (width > 0 && width + length + 2*writer->getIndent() > 60) {
			width = 0;
			*writer << endl;
			if (!extraIndent) {
				writer->indent();
				extraIndent = true;
			}
		} else {
			width += length;
			*writer << " ";
		}
		*writer << key.c_str() << "=\"" << value.c_str() << "\"";
	}
	attributes.clear();
	if (extraIndent) writer->outdent();
}

string XMLWriter::normalize(string s) {
    string str = "";
    char buffer[20];

    int len = s.length();
    for (int i = 0; i < len; i++) {
        char ch = s[i];
        switch (ch) {
            case '<': {
                str.append("&lt;");
                break;
            }
            case '>': {
                str.append("&gt;");
                break;
            }
            case '&': {
                str.append("&amp;");
                break;
            }
            case '"': {
                str.append("&quot;");
                break;
            }
            case '\r':
            case '\n': {
                sprintf(buffer, "&#%ud", ch);
                str.append(buffer);
                str.append(";");
                break;
            }
            default: {
//                if (ch > 0x00FF) {
//                    sprintf(buffer, "&#x%4.4x", ch);
//                    str.append(buffer);
//                    str.append(";");
//                } else {
                    str.append(&ch, 1);
//                }
            }
        }
    }

    return str;
}

string XMLWriter::normalizeText(string s) {
    string str = "";

    int len = s.length();
    for (int i = 0; i < len; i++) {
        char ch = s[i];
        switch (ch) {
            case '<': {
                str.append("&lt;");
                break;
            }
            case '>': {
                str.append("&gt;");
                break;
            }
            case '&': {
                str.append("&amp;");
                break;
            }
            default: {
//                if (ch > 0x00FF) {
//                    sprintf(buffer, "&#x%4.4x", ch);
//                    str.append(buffer);
//                    str.append(";");
//                } else {
                    str.append(&ch, 1);
//                }
            }
        }
    }
    return str;
}

void XMLWriter::checkNameValid(string) {
// Could be added.
//    if (!XMLCharacterProperties.validName(s)) throw new RuntimeException("Invalid name: "+s);
}


} // cheprep
