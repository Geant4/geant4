// Copyright FreeHEP, 2005.

#include <iostream>

#include "cheprep/DefaultHepRepTypeTree.h"

using namespace std;
using namespace HEPREP;

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepTypeTree.cc,v 1.7 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

DefaultHepRepTypeTree::DefaultHepRepTypeTree(HepRepTreeID* typeTree)
    : DefaultHepRepTreeID(typeTree->getName(), typeTree->getVersion()) {
    delete typeTree;
}

DefaultHepRepTypeTree::~DefaultHepRepTypeTree() {
    for (vector<HepRepType*>::iterator i1 = types.begin(); i1 != types.end(); i1++) {
        delete (*i1);
    }
}

HepRepTypeTree* DefaultHepRepTypeTree::copy() {
    cerr << "DefaultHepRepTypeTree::copy() not implemented." << endl;
    return NULL;
}

void DefaultHepRepTypeTree::addType(HepRepType* type) {
    // FIXME should check if type already exists
    types.push_back(type);
}

vector<HepRepType*> DefaultHepRepTypeTree::getTypeList() {
    return types;
}

HepRepType* DefaultHepRepTypeTree::getType(string /*typeName*/) {
    cerr << "DefaultHepRepTypeTree::getType(string) not implemented." << endl;
    return NULL;
}

} // cheprep
