
#include "DefaultHepRepInstanceTree.h"

using namespace std;
using namespace HEPREP;

DefaultHepRepInstanceTree::DefaultHepRepInstanceTree(string name, string version, HepRepTreeID* typeTree)
    : DefaultHepRepTreeID(name, version), typeTree(typeTree) {
}

DefaultHepRepInstanceTree::~DefaultHepRepInstanceTree() {
    for (vector<HepRepInstance*>::iterator i1 = instances.begin(); i1 != instances.end(); i1++) {
        delete (*i1);
    }
    instances.clear();
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

