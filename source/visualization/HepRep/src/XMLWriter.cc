//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
#include "XMLWriter.h"

using namespace std;

XMLWriter::XMLWriter(ostream* out, string indentString, string defaultNameSpace)
    : defaultNameSpace(defaultNameSpace) {
    writer = new IndentPrintWriter(out);
    writer->setIndentString("  ");
    closed = true;
    dtdName = NULL;
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

    closed = false;
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
    if (dtdName != NULL) {
        cerr << "XMLWriter::ReferToDTD cannot be called twice" << endl;
    }
    dtdName = &name;
    *writer << "<!DOCTYPE " << name.c_str() << " PUBLIC \"" << pid.c_str() << "\" \"" << ref.c_str() << "\">" << endl;
}

void XMLWriter::referToDTD(string name, string system) {
    if (dtdName != NULL) {
        cerr << "XMLWriter::ReferToDTD cannot be called twice";
    }
    dtdName = &name;
    *writer << "<!DOCTYPE " << name.c_str() << " SYSTEM \"" << system.c_str() << "\">" << endl;
}

void XMLWriter::closeDoc() {
    if (!closed) {
        if (!openTags.empty()) {
            cerr << "Not all tags were closed before closing XML document:" << endl;
            while (!openTags.empty()) {
                cerr << "   </" << openTags.top().c_str() << ">" << endl;
                openTags.pop();
            }
        }
        closed = true;
    }
}

void XMLWriter::printComment(string comment) {
    if (comment.find("--") >= 0) {
        cerr << "XMLWriter::printComment '--' sequence not allowed in comment" << endl;
    }
    *writer << "<!--" << normalizeText(comment).c_str() << "-->" << endl;
}

void XMLWriter::print(string text) {
    *writer << normalizeText(text).c_str();
}

void XMLWriter::println(string text) {
    print(text);
    *writer << endl;
}

void XMLWriter::openTag(string ns, string name) {
    if (ns.compare(defaultNameSpace) == 0) {
        openTag(name);
    } else {
        openTag(ns.append(":").append(name));
    }
}

void XMLWriter::openTag(string name) {
    checkNameValid(name);
    if (openTags.empty() && dtdName != NULL && dtdName->compare(name)) {
        cerr << "XMLWriter::openTag(), First tag: '" << name.c_str() << "' not equal to DTD id: '" << dtdName->c_str() << "'" << endl;
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

void XMLWriter::printTag(string ns, string name) {
    if (ns.compare(defaultNameSpace) == 0) {
        printTag(name);
    } else {
        printTag(ns.append(":").append(name));
    }
}

void XMLWriter::printTag(string name) {
    checkNameValid(name);
    *writer << "<" << name.c_str();
    printAttributes(name.length());
    *writer << "/>" << endl;
}

void XMLWriter::setAttribute(string name, string value) {
    attributes[name] = value;
}

void XMLWriter::setAttribute(string ns, string name, string value) {
    attributes[ns.append(":").append(name)] = value;
}

void XMLWriter::setAttribute(string name, double value) {
    char buffer[40];
    sprintf(buffer, "%f", value);
    attributes[name] = buffer;
}

void XMLWriter::setAttribute(string ns, string name, double value) {
    char buffer[40];
    sprintf(buffer, "%f", value);
    attributes[ns.append(":").append(name)] = buffer;
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

void XMLWriter::checkNameValid(string s) {
// Could be added.
//    if (!XMLCharacterProperties.validName(s)) throw new RuntimeException("Invalid name: "+s);
}
