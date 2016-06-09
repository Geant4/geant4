
#include <iostream>

#include "DefaultHepRepInstance.h"

using namespace std;
using namespace HEPREP;

DefaultHepRepInstance::DefaultHepRepInstance(HepRepInstance* instance, HepRepType* heprepType)
    : DefaultHepRepAttribute(), parent(instance), type(heprepType) {

    if (type == NULL) cerr << "HepRepInstance cannot be created without a HepRepType." << endl;
    ((HepRepInstance*)parent)->addInstance(this);
}

DefaultHepRepInstance::DefaultHepRepInstance(HepRepInstanceTree* instanceTree, HepRepType* heprepType)
    : DefaultHepRepAttribute(), parent(instanceTree), type(heprepType) {

    if (type == NULL) cerr << "HepRepInstance cannot be created without a HepRepType." << endl;
    ((HepRepInstanceTree*)parent)->addInstance(this);
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

HepRepInstance* DefaultHepRepInstance::copy(HepRep*, HepRepInstance*, HepRepSelectFilter*) {
    cerr << "DefaultHepRepInstance::copy(HepRep*, HepRepInstance*, HepRepSelectFilter*) not implemented." << endl;
    return NULL;
}

HepRepInstance* DefaultHepRepInstance::copy(HepRep*, HepRepInstanceTree*, HepRepSelectFilter*) {
    cerr << "DefaultHepRepInstance::copy(HepRep*, HepRepInstanceTree*, HepRepSelectFilter*) not implemented." << endl;
    return NULL;
}

HepRepType* DefaultHepRepInstance::getType() {
    return type;
}

bool DefaultHepRepInstance::addPoint(HepRepPoint* point) {
    points.push_back(point);
    return true;
}

vector<HepRepPoint*>* DefaultHepRepInstance::getPoints() {
    return &points;
}

bool DefaultHepRepInstance::addInstance(HepRepInstance* instance) {
    instances.push_back(instance);
    return true;
}

void DefaultHepRepInstance::removeInstance(HepRepInstance*) {
    cerr << "DefaultHepRepInstance::removeInstance(HepRepInstance*) not implemented." << endl;
}

vector<HepRepInstance*>* DefaultHepRepInstance::getInstances() {
    return &instances;
}

HepRepAttValue* DefaultHepRepInstance::getAttValue(string name) {
    HepRepAttValue* value = getAttValueFromNode(name);
    return (value != NULL) ? value : type->getAttValue(name);
}
