
#include <iostream>

#include "DefaultHepRepTypeTree.h"

using namespace std;
using namespace HEPREP;

DefaultHepRepTypeTree::DefaultHepRepTypeTree(HepRepTreeID* typeTree)
    : DefaultHepRepTreeID(typeTree->getName(), typeTree->getVersion()) {
    delete typeTree;
}

DefaultHepRepTypeTree::~DefaultHepRepTypeTree() {
    for (set<HepRepType*>::iterator i1 = types.begin(); i1 != types.end(); i1++) {
        delete (*i1);
    }
}

HepRepTypeTree* DefaultHepRepTypeTree::copy() {
    cerr << "DefaultHepRepTypeTree::copy() not implemented." << endl;
    return NULL;
}

void DefaultHepRepTypeTree::addType(HepRepType* type) {
    types.insert(type);
}

set<HepRepType*> DefaultHepRepTypeTree::getTypes() {
    return types;
}

HepRepType* DefaultHepRepTypeTree::getType(string /*typeName*/) {
    cerr << "DefaultHepRepTypeTree::getType(string) not implemented." << endl;
    return NULL;
}
