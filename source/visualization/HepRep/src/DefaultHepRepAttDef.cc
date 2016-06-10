// Copyright FreeHEP, 2005.

#include <iostream>
#include <cstring>
#include <cctype>
#include <algorithm>

#include "cheprep/DefaultHepRepAttDef.h"

using namespace std;
using namespace HEPREP;

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepAttDef.cc 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

DefaultHepRepAttDef::DefaultHepRepAttDef(string aName, string aDesc, string aCategory, string anExtra)
    : name(aName), desc(aDesc), category(aCategory), extra(anExtra) {
}

DefaultHepRepAttDef::~DefaultHepRepAttDef() {
}

HepRepAttDef* DefaultHepRepAttDef::copy() {
    return new DefaultHepRepAttDef(name, desc, category, extra);
}

string DefaultHepRepAttDef::getName() {
    return name;
}

string DefaultHepRepAttDef::getLowerCaseName() {
    string s = name;
    transform(s.begin(), s.end(), s.begin(), (int(*)(int)) tolower);
    return s;
}

string DefaultHepRepAttDef::getDescription() {
    return desc;
}

string DefaultHepRepAttDef::getCategory() {
    return category;
}

string DefaultHepRepAttDef::getExtra() {
    return extra;
}

} // cheprep
