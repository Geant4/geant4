// Copyright FreeHEP, 2005.

#include <iostream>

#include "cheprep/DefaultHepRepInstance.h"

using namespace std;
using namespace HEPREP;

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepInstance.cc,v 1.7 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

DefaultHepRepInstance::DefaultHepRepInstance(HepRepInstance* instance, HepRepType* heprepType)
    : DefaultHepRepAttribute(), parent(instance), type(heprepType) {

    if (type == NULL) cerr << "HepRepInstance cannot be created without a HepRepType." << endl;
    parent->addInstance(this);
}

DefaultHepRepInstance::DefaultHepRepInstance(HepRepInstanceTree* instanceTree, HepRepType* heprepType)
    : DefaultHepRepAttribute(), parent(NULL), type(heprepType) {

    if (type == NULL) cerr << "HepRepInstance cannot be created without a HepRepType." << endl;
    instanceTree->addInstance(this);
}

DefaultHepRepInstance::~DefaultHepRepInstance() {
    parent = NULL;
    type = NULL;
    for (vector<HepRepInstance*>::iterator i1 = instances.begin(); i1 != instances.end(); i1++) {
        delete (*i1);
    }
    for (vector<HepRepPoint*>::iterator i2 = points.begin(); i2 != points.end(); i2++) {
        delete (*i2);
    }
}

void DefaultHepRepInstance::overlay(HepRepInstance *) {
    cerr << "DefaultHepRepInstance::overlay(HepRepInstance * instance) not implemented." << endl;
}

HepRepInstance* DefaultHepRepInstance::copy(HepRepTypeTree*, HepRepInstance*, HepRepSelectFilter*) {
    cerr << "DefaultHepRepInstance::copy(HepRepTypeTree*, HepRepInstance*, HepRepSelectFilter*) not implemented." << endl;
    return NULL;
}

HepRepInstance* DefaultHepRepInstance::copy(HepRepTypeTree*, HepRepInstanceTree*, HepRepSelectFilter*) {
    cerr << "DefaultHepRepInstance::copy(HepRepTypeTree*, HepRepInstanceTree*, HepRepSelectFilter*) not implemented." << endl;
    return NULL;
}

HepRepType* DefaultHepRepInstance::getType() {
    return type;
}

void DefaultHepRepInstance::addPoint(HepRepPoint* point) {
    points.push_back(point);
}

vector<HepRepPoint*> DefaultHepRepInstance::getPoints() {
    return points;
}

HepRepInstance* DefaultHepRepInstance::getSuperInstance() {
    return parent;
}

void DefaultHepRepInstance::addInstance(HepRepInstance* instance) {
    instances.push_back(instance);
}

void DefaultHepRepInstance::removeInstance(HepRepInstance*) {
    cerr << "DefaultHepRepInstance::removeInstance(HepRepInstance*) not implemented." << endl;
}

vector<HepRepInstance*> DefaultHepRepInstance::getInstances() {
    return instances;
}

HepRepAttValue* DefaultHepRepInstance::getAttValue(string name) {
    HepRepAttValue* value = getAttValueFromNode(name);
    return (value != NULL) ? value : type->getAttValue(name);
}

} // cheprep
