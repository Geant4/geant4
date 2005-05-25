// Copyright FreeHEP, 2005.

#include "cheprep/DefaultHepRepInstanceTree.h"

using namespace std;
using namespace HEPREP;

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRepInstanceTree.cc,v 1.10 2005-05-25 23:22:25 duns Exp $
 */
namespace cheprep {

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

HepRepInstanceTree* DefaultHepRepInstanceTree::copy(HepRepTypeTree*, HepRepSelectFilter*) {
    cerr << "DefaultHepRepInstanceTree::copy(HepRepTypeTree*, HepRepSelectFilter*) not implemented." << endl;
    return NULL;
}

void DefaultHepRepInstanceTree::addInstance(HepRepInstance* instance) {
    instances.push_back(instance);
}

void DefaultHepRepInstanceTree::removeInstance(HepRepInstance*) {
    cerr << "DefaultHepRepInstanceTree::removeInstance(HepRepInstance*) not implemented." << endl;
}

vector<HepRepInstance*> DefaultHepRepInstanceTree::getInstances() {
    return instances;
}

void DefaultHepRepInstanceTree::addInstanceTree(HepRepTreeID* treeID) {
    instanceTrees.push_back(treeID);
}

HepRepTreeID* DefaultHepRepInstanceTree::getTypeTree() {
    return typeTree;
}

vector<HepRepTreeID*> DefaultHepRepInstanceTree::getInstanceTreeList() {
    return instanceTrees;
}

} // cheprep

