
#include "DefaultHepRepInstanceTree.h"

using namespace std;
using namespace HEPREP;

DefaultHepRepInstanceTree::DefaultHepRepInstanceTree(string name, string version, HepRepTreeID* typeTree)
    : DefaultHepRepTreeID(name, version), typeTree(typeTree) {
}

DefaultHepRepInstanceTree::~DefaultHepRepInstanceTree() {
cout << "Delete Tree " << getName() << endl;
    for (vector<HepRepInstance*>::iterator i1 = instances.begin(); i1 != instances.end(); i1++) {
        delete (*i1);
    }
    instances.clear();

    for (vector<HepRepTreeID*>::iterator i2 = instanceTrees.begin(); i2 != instanceTrees.end(); i2++) {
        delete (*i2);
    }
    instanceTrees.clear();
}

void DefaultHepRepInstanceTree::overlay(HepRepInstanceTree *) {
    cerr << "DefaultHepRepInstanceTree::overlay(HepRepInstanceTree * instanceTree) not implemented." << endl;
}

HepRepTreeID* DefaultHepRepInstanceTree::copy() {
    return DefaultHepRepTreeID::copy();
}

HepRepInstanceTree* DefaultHepRepInstanceTree::copy(HepRep*, HepRepSelectFilter*) {
    cerr << "DefaultHepRepInstanceTree::copy(HepRep*, HepRepSelectFilter*) not implemented." << endl;
    return NULL;
}

bool DefaultHepRepInstanceTree::addInstance(HepRepInstance* instance) {
    instances.push_back(instance);
    return true;
}

void DefaultHepRepInstanceTree::removeInstance(HepRepInstance*) {
    cerr << "DefaultHepRepInstanceTree::removeInstance(HepRepInstance*) not implemented." << endl;
}

vector<HepRepInstance*>* DefaultHepRepInstanceTree::getInstances() {
    return &instances;
}

bool DefaultHepRepInstanceTree::addInstanceTree(HepRepTreeID* treeID) {
    instanceTrees.push_back(treeID);
    return true;
}

HepRepTreeID* DefaultHepRepInstanceTree::getTypeTree() {
    return typeTree;
}

vector<HepRepTreeID*>* DefaultHepRepInstanceTree::getInstanceTrees() {
    return &instanceTrees;
}

