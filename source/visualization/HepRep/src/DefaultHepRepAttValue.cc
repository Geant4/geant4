
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cstdio>
#include <iostream>

#include "HEPREP/HepRepConstants.h"

#include "DefaultHepRepAttValue.h"

using namespace std;
using namespace HEPREP;

DefaultHepRepAttValue::DefaultHepRepAttValue(string name, string value, int showLabel)
    : name(name), type(HepRepConstants::TYPE_STRING), stringValue(value), showLabelValue(showLabel) {

    init();
}

DefaultHepRepAttValue::DefaultHepRepAttValue(string name, long value, int showLabel)
    : name(name), type(HepRepConstants::TYPE_LONG), longValue(value), showLabelValue(showLabel) {

    init();
}

DefaultHepRepAttValue::DefaultHepRepAttValue(string name, int value, int showLabel)
    : name(name), type(HepRepConstants::TYPE_INT), longValue(value), showLabelValue(showLabel) {

    init();
}

DefaultHepRepAttValue::DefaultHepRepAttValue(string name, double value, int showLabel)
    : name(name), type(HepRepConstants::TYPE_DOUBLE), doubleValue(value), showLabelValue(showLabel) {

    init();
}

DefaultHepRepAttValue::DefaultHepRepAttValue(string name, bool value, int showLabel)
    : name(name), type(HepRepConstants::TYPE_BOOLEAN), booleanValue(value), showLabelValue(showLabel) {

    init();
}

DefaultHepRepAttValue::DefaultHepRepAttValue(string name, vector<double> value, int showLabel)
    : name(name), type(HepRepConstants::TYPE_COLOR), colorValue(value), showLabelValue(showLabel) {

    init();
}

DefaultHepRepAttValue::~DefaultHepRepAttValue() {
}

void DefaultHepRepAttValue::init() {
    labelStrings[0] = "NAME";
    labelStrings[1] = "DESC";
    labelStrings[2] = "VALUE";
    labelStrings[3] = "EXTRA";
}

HepRepAttValue* DefaultHepRepAttValue::copy() {
    return NULL;
}

string DefaultHepRepAttValue::getName() {
    return name;
}

string DefaultHepRepAttValue::getLowerCaseName() {
    char *tmp = new char[strlen(name.c_str())];
    strcpy(tmp, name.c_str());
    int i = -1;
    do {
        i++;
        tmp[i] = tolower(tmp[i]);
    } while (tmp[i] != 0);
    return tmp;
}

int DefaultHepRepAttValue::getType() {
    return type;
}

string DefaultHepRepAttValue::getTypeName() {
    switch(type) {
        case HepRepConstants::TYPE_COLOR: return("Color");
        case HepRepConstants::TYPE_STRING: return("String");
        case HepRepConstants::TYPE_LONG: return("long");
        case HepRepConstants::TYPE_INT: return("int");
        case HepRepConstants::TYPE_DOUBLE: return("double");
        case HepRepConstants::TYPE_BOOLEAN: return("boolean");
        default: return "Unknown type stored in HepRepAttDef";
    }
}

int DefaultHepRepAttValue::showLabel() {
    return showLabelValue;
}

string DefaultHepRepAttValue::getString() {
    if (type != HepRepConstants::TYPE_STRING) cerr << "Trying to access AttValue as 'string'" << endl;
    return stringValue;
}

long DefaultHepRepAttValue::getLong() {
    if (type != HepRepConstants::TYPE_LONG) cerr << "Trying to access AttValue as 'long'" << endl;
    return longValue;
}

int DefaultHepRepAttValue::getInteger() {
    if (type != HepRepConstants::TYPE_INT) cerr << "Trying to access AttValue as 'int'" << endl;
    return (int)longValue;
}

double DefaultHepRepAttValue::getDouble() {
    if (type != HepRepConstants::TYPE_DOUBLE) cerr << "Trying to access AttValue as 'double'" << endl;
    return doubleValue;
}

bool DefaultHepRepAttValue::getBoolean() {
    if (type != HepRepConstants::TYPE_BOOLEAN) cerr << "Trying to access AttValue as 'boolean'" << endl;
    return booleanValue;
}

vector<double> DefaultHepRepAttValue::getColor() {
    if (type != HepRepConstants::TYPE_COLOR) cerr << "Trying to access AttValue as 'color'" << endl;
    return colorValue;
}


string DefaultHepRepAttValue::getAsString() {
    char buffer[40];
    switch(type) {
        case HepRepConstants::TYPE_COLOR:   sprintf(buffer, "%4.2f, %4.2f, %4.2f, %4.2f",
                                                colorValue[0],
                                                colorValue[1],
                                                colorValue[2],
                                                (colorValue.size() > 3) ? colorValue[3] : 1.0);
                                            return buffer;
        case HepRepConstants::TYPE_STRING:  return getString();
        case HepRepConstants::TYPE_LONG:    sprintf(buffer, "%ld", getLong());
                                            return buffer;
        case HepRepConstants::TYPE_INT:     sprintf(buffer, "%d", getInteger());
                                            return buffer;
        case HepRepConstants::TYPE_DOUBLE:  sprintf(buffer, "%f", getDouble());
                                            return buffer;
        case HepRepConstants::TYPE_BOOLEAN: return (getBoolean()) ? "true" : "false";
        default:                            return "Unknown typecode";
    }
}



string DefaultHepRepAttValue::toShowLabel() {
    string label = "";
    bool first = true;
    if (showLabel() == HepRepConstants::SHOW_NONE) {
        label = "NONE";
    } else {
        for (int i=0; i<16; i++) {
            if (((showLabel() >> i) & 0x0001) == 0x0001) {
                if (first) {
                    first = false;
                } else {
                    label.append(", ");
                }
                if (i < LABELSTRINGS_LEN) {
                    label.append(labelStrings[i]);
                } else {
                    char hex[20];
                    sprintf(hex, "%0x", 1 << i);
                    label.append(hex);
                }
            }
        }
    }
    return label;
}
