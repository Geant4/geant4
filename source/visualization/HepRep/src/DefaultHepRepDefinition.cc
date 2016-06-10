// Copyright FreeHEP, 2005.

#include "cheprep/DefaultHepRepDefinition.h"
#include "cheprep/DefaultHepRepAttDef.h"

#include <iostream>
#include <algorithm>

using namespace std;
using namespace HEPREP;

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepDefinition.cc 66373 2012-12-18 09:41:34Z gcosmo $
 */
namespace cheprep {

DefaultHepRepDefinition::DefaultHepRepDefinition()
    : DefaultHepRepAttribute() {
}

DefaultHepRepDefinition::~DefaultHepRepDefinition() {
    set<HepRepAttDef *> list = getAttDefsFromNode();
    for (set<HepRepAttDef*>::iterator i1 = list.begin(); i1 != list.end(); i1++) {
        delete (*i1);
    }
}

set<HepRepAttDef *> DefaultHepRepDefinition::getAttDefsFromNode() {
    set<HepRepAttDef*> attSet;
    for (map<string, HepRepAttDef*>::iterator i = attDefs.begin(); i != attDefs.end(); i++) {
        attSet.insert((*i).second);
    }
    return attSet;
}

void DefaultHepRepDefinition::addAttDef(HepRepAttDef* hepRepAttDef) {
    string lowerCaseName = hepRepAttDef->getLowerCaseName();
    if (attDefs[lowerCaseName] != NULL) delete attDefs[lowerCaseName];
    attDefs[lowerCaseName] = hepRepAttDef;
}

void DefaultHepRepDefinition::addAttDef(string name, string desc, string type, string extra) {
    addAttDef(new DefaultHepRepAttDef(name, desc, type, extra));
}

HepRepAttDef* DefaultHepRepDefinition::getAttDefFromNode(string name) {
    string s = name;
    transform(s.begin(), s.end(), s.begin(), (int(*)(int)) tolower);
    return (attDefs.count(s) > 0) ? attDefs[s] : NULL;    
}

} // cheprep
