
#include "DefaultHepRepDefinition.h"
#include "DefaultHepRepAttDef.h"

#include <iostream>

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
    attDefs[hepRepAttDef->getLowerCaseName()] = hepRepAttDef;
    return true;
}

bool DefaultHepRepDefinition::addAttDef(string name, string desc, string type, string extra) {
    return addAttDef(new DefaultHepRepAttDef(name, desc, type, extra));
}

HepRepAttDef* DefaultHepRepDefinition::getAttDefFromNode(string name) {
    map<string, HepRepAttDef*>::iterator i = attDefs.find(name);
    return i != attDefs.end() ? (*i).second : NULL;
}

