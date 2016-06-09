
#include "DefaultHepRepDefinition.h"
#include "DefaultHepRepAttDef.h"

#include <iostream>
#include <algorithm>

using namespace std;
using namespace HEPREP;

DefaultHepRepDefinition::DefaultHepRepDefinition()
    : DefaultHepRepAttribute() {
}

DefaultHepRepDefinition::~DefaultHepRepDefinition() {
    vector<HepRepAttDef *>* list = getAttDefsFromNode();
    for (vector<HepRepAttDef*>::iterator i1 = list->begin(); i1 != list->end(); i1++) {
        delete (*i1);
    }
}

vector<HepRepAttDef *>* DefaultHepRepDefinition::getAttDefsFromNode() {
    attList.clear();
    for (map<string, HepRepAttDef*>::iterator i = attDefs.begin(); i != attDefs.end(); i++) {
        attList.push_back((*i).second);
    }
    return &attList;
}

bool DefaultHepRepDefinition::addAttDef(HepRepAttDef* hepRepAttDef) {
    string lowerCaseName = hepRepAttDef->getLowerCaseName();
    if (attDefs[lowerCaseName] != NULL) delete attDefs[lowerCaseName];
    attDefs[lowerCaseName] = hepRepAttDef;
    return true;
}

bool DefaultHepRepDefinition::addAttDef(string name, string desc, string type, string extra) {
    return addAttDef(new DefaultHepRepAttDef(name, desc, type, extra));
}

HepRepAttDef* DefaultHepRepDefinition::getAttDefFromNode(string name) {
    string s = name;
    transform(s.begin(), s.end(), s.begin(), (int(*)(int)) tolower);
    return (attDefs.count(s) > 0) ? attDefs[s] : NULL;    
}

