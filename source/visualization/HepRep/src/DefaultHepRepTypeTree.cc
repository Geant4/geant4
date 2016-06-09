
#include <iostream>

#include "DefaultHepRepTypeTree.h"

using namespace std;
using namespace HEPREP;

DefaultHepRepTypeTree::DefaultHepRepTypeTree(HepRepTreeID* typeTree)
    : DefaultHepRepTreeID(typeTree->getName(), typeTree->getVersion()) {
    delete typeTree;
}

DefaultHepRepTypeTree::~DefaultHepRepTypeTree() {
    for (vector<HepRepType*>::iterator i1 = types.begin(); i1 != types.end(); i1++) {
        delete (*i1);
    }
}

HepRepTreeID* DefaultHepRepTypeTree::copy() {
    return DefaultHepRepTreeID::copy();
}

HepRepTypeTree* DefaultHepRepTypeTree::copy(HepRep*) {
    cerr << "DefaultHepRepTypeTree::copy(HepRep*) not implemented." << endl;
    return NULL;
}

bool DefaultHepRepTypeTree::addType(HepRepType* type) {
    types.push_back(type);
    return true;
}

vector<HepRepType*>* DefaultHepRepTypeTree::getTypes() {
    return &types;
}

