
#include "DefaultHepRep.h"

using namespace std;
using namespace HEPREP;

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

vector<string>* DefaultHepRep::getLayerOrder() {
    return &layers;
}

bool DefaultHepRep::addLayer(string layer) {
    layers.push_back(layer);
    return true;
}

bool DefaultHepRep::addTypeTree(HepRepTypeTree* typeTree) {
    typeTrees.push_back(typeTree);
    return true;
}

void DefaultHepRep::removeTypeTree(HepRepTypeTree*) {
    cerr << "DefaultHepRep::removeTypeTree(HepRepTypeTree*) not implemented." << endl;
}

HepRepTypeTree* DefaultHepRep::getTypeTree(string, string) {
    cerr << "DefaultHepRep::getTypeTree(string, string) not implemented." << endl;
    return NULL;
}

vector<HepRepTypeTree*>* DefaultHepRep::getTypeTrees() {
    return &typeTrees;
}

HepRepType* DefaultHepRep::getType(string) {
    cerr << "DefaultHepRep::getType(string) not implemented." << endl;
    return NULL;
}

bool DefaultHepRep::addInstanceTree(HepRepInstanceTree* instanceTree) {
    instanceTrees.push_back(instanceTree);
    return true;
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

vector<HepRepInstanceTree*>* DefaultHepRep::getInstanceTrees() {
    return &instanceTrees;
}

