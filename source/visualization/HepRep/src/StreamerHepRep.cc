
#include "StreamerHepRep.h"

using namespace std;
using namespace HEPREP;

StreamerHepRep::StreamerHepRep(HepRepWriter* streamer) {
    streamer->write(this);
}

StreamerHepRep::~StreamerHepRep() {
}

HepRep* StreamerHepRep::copy(HepRepSelectFilter*) {
    return NULL;
}

vector<string>* StreamerHepRep::getLayerOrder() {
    return &layers;
}

bool StreamerHepRep::addLayer(string layer) {
    layers.push_back(layer);
    return true;
}

bool StreamerHepRep::addTypeTree(HepRepTypeTree*) {
    return true;
}

void StreamerHepRep::removeTypeTree(HepRepTypeTree*) {
}

HepRepTypeTree* StreamerHepRep::getTypeTree(string, string) {
    return NULL;
}

vector<HepRepTypeTree*>* StreamerHepRep::getTypeTrees() {
    return NULL;
}

HepRepType* StreamerHepRep::getType(string) {
    return NULL;
}

bool StreamerHepRep::addInstanceTree(HepRepInstanceTree*) {
    return true;
}

void StreamerHepRep::removeInstanceTree(HepRepInstanceTree*) {
}

HepRepInstanceTree* StreamerHepRep::getInstanceTreeTop(string, string) {
    return NULL;
}

HepRepInstanceTree* StreamerHepRep::getInstances(string, string,
                                       vector<string>) {
    return NULL;
}

HepRepInstanceTree* StreamerHepRep::getInstancesAfterAction(
                                string,
                                string,
                                vector<string>,
                                vector<HepRepAction*>,
                                bool,
                                bool,
                                bool,
                                vector<string>) {
    return NULL;
}

string StreamerHepRep::checkForException() {
    return NULL;
}

vector<HepRepInstanceTree*>* StreamerHepRep::getInstanceTrees() {
    return NULL;
}

