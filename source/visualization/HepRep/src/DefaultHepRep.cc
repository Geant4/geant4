// Copyright FreeHEP, 2005.

#include "cheprep/DefaultHepRep.h"

using namespace std;
using namespace HEPREP;

/**
 * @author Mark Donszelmann
 * @version $Id: DefaultHepRep.cc,v 1.9 2005-06-02 21:28:45 duns Exp $
 */
namespace cheprep {

DefaultHepRep::DefaultHepRep() {
}

DefaultHepRep::~DefaultHepRep() {
    for (vector<HepRepTypeTree*>::iterator i1 = typeTrees.begin(); i1 != typeTrees.end(); i1++) {
        delete (*i1);
    }
    for (vector<HepRepInstanceTree*>::iterator i2 = instanceTrees.begin(); i2 != instanceTrees.end(); i2++) {
        delete (*i2);
    }
}

HepRep* DefaultHepRep::copy(HepRepSelectFilter*) {
    cerr << "DefaultHepRep::copy(HepRepSelectFilter*) not implemented." << endl;
    return NULL;
}

vector<string> DefaultHepRep::getLayerOrder() {
    return layers;
}

void DefaultHepRep::addLayer(string layer) {
    layers.push_back(layer);
}

void DefaultHepRep::addTypeTree(HepRepTypeTree* typeTree) {
    // FIXME check if already exist    
    typeTrees.push_back(typeTree);
}

void DefaultHepRep::removeTypeTree(HepRepTypeTree*) {
    cerr << "DefaultHepRep::removeTypeTree(HepRepTypeTree*) not implemented." << endl;
}

HepRepTypeTree* DefaultHepRep::getTypeTree(string, string) {
    cerr << "DefaultHepRep::getTypeTree(string, string) not implemented." << endl;
    return NULL;
}

vector<HepRepTypeTree*> DefaultHepRep::getTypeTreeList() {
    return typeTrees;
}

void DefaultHepRep::addInstanceTree(HepRepInstanceTree* instanceTree) {
    // FIXME check if already exist    
    instanceTrees.push_back(instanceTree);
}

void DefaultHepRep::overlayInstanceTree(HepRepInstanceTree *) {
    cerr << "DefaultHepRep::overlayInstanceTree(HepRepInstanceTree * instanceTree) not implemented." << endl;
}

void DefaultHepRep::removeInstanceTree(HepRepInstanceTree*) {
    cerr << "DefaultHepRep::removeInstanceTree(HepRepInstanceTree*) not implemented." << endl;
}

HepRepInstanceTree* DefaultHepRep::getInstanceTreeTop(string, string) {
    cerr << "DefaultHepRep::getInstanceTreeTop(string, string) not implemented." << endl;
    return NULL;
}

HepRepInstanceTree* DefaultHepRep::getInstances(string, string,
                                                vector<string>) {
    cerr << "DefaultHepRep::getInstances(string, string, vector<string>) not implemented." << endl;
    return NULL;
}

HepRepInstanceTree* DefaultHepRep::getInstancesAfterAction(
                                string,
                                string,
                                vector<string>,
                                vector<HepRepAction*>,
                                bool,
                                bool,
                                bool,
                                vector<string>) {
    cerr << "DefaultHepRep::getInstancesAfterAction(string, string, vector<string>, vector<HepRepAction*>, bool, bool, bool, vector<string>) not implemented." << endl;
    return NULL;
}

string DefaultHepRep::checkForException() {
    return NULL;
}

vector<HepRepInstanceTree*> DefaultHepRep::getInstanceTreeList() {
    return instanceTrees;
}

} // cheprep

