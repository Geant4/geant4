
#include <cstring>
#include <cctype>

#include "DefaultHepRepAttDef.h"

using namespace std;
using namespace HEPREP;

DefaultHepRepAttDef::DefaultHepRepAttDef(string name, string desc, string category, string extra)
    : name(name), desc(desc), category(category), extra(extra) {
}

DefaultHepRepAttDef::~DefaultHepRepAttDef() {
}

HepRepAttDef* DefaultHepRepAttDef::copy() {
    return NULL;
}

string DefaultHepRepAttDef::getName() {
    return name;
}

string DefaultHepRepAttDef::getLowerCaseName() {
    char* tmp = new char[strlen(name.c_str())];
    strcpy(tmp, name.c_str());
    int i = -1;
    do {
        i++;
        tmp[i] = tolower(tmp[i]);
    } while (tmp[i] != 0);
    return tmp;
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

