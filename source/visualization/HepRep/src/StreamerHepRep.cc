
#include "StreamerHepRep.h"

using namespace std;
using namespace HEPREP;

StreamerHepRep::StreamerHepRep(HepRepWriter* streamer) {
    streamer->write(this);
}

StreamerHepRep::~StreamerHepRep() {
}

HepRep* StreamerHepRep::copy(HepRepSelectFilter* filter) {
    return NULL;
}

vector<string>* StreamerHepRep::getLayerOrder() {
    return &layers;
}

bool StreamerHepRep::addLayer(string layer) {
    layers.push_back(layer);
    return true;
}

bool StreamerHepRep::addTypeTree(HepRepTypeTree* typeTree) {
    return true;
}

void StreamerHepRep::removeTypeTree(HepRepTypeTree* typeTree) {
}

HepRepTypeTree* StreamerHepRep::getTypeTree(string name, string version) {
    return NULL;
}

vector<HepRepTypeTree*>* StreamerHepRep::getTypeTrees() {
    return NULL;
}

HepRepType* StreamerHepRep::getType(string name) {
    return NULL;
}

bool StreamerHepRep::addInstanceTree(HepRepInstanceTree* instanceTree) {
    return true;
}

void StreamerHepRep::removeInstanceTree(HepRepInstanceTree* instanceTree) {
}

HepRepInstanceTree* StreamerHepRep::getInstanceTreeTop(string name, string version) {
    return NULL;
}

HepRepInstanceTree* StreamerHepRep::getInstances(string name, string version,
                                       vector<string> typeNames) {
    return NULL;
}

HepRepInstanceTree* StreamerHepRep::getInstancesAfterAction(
                                string name,
                                string version,
                                vector<string> typeNames,
                                vector<HepRepAction*> actions,
                                bool getPoints,
                                bool getDrawAtts,
                                bool getNonDrawAtts,
                                vector<string> invertAtts) {
    return NULL;
}

string StreamerHepRep::checkForException() {
    return NULL;
}

vector<HepRepInstanceTree*>* StreamerHepRep::getInstanceTrees() {
    return NULL;
}

