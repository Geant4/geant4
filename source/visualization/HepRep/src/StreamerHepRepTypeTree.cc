
#include <iostream>

#include "StreamerHepRepTypeTree.h"

using namespace std;
using namespace HEPREP;

StreamerHepRepTypeTree::StreamerHepRepTypeTree(HepRepWriter* streamer, HepRepTreeID* typeTree)
    : DefaultHepRepTreeID(typeTree->getName(), typeTree->getVersion()) {

    streamer->write(this);
}

StreamerHepRepTypeTree::~StreamerHepRepTypeTree() {
}

HepRepTreeID* StreamerHepRepTypeTree::copy() {
    return DefaultHepRepTreeID::copy();
}

HepRepTypeTree* StreamerHepRepTypeTree::copy(HepRep* heprep) {
    return NULL;
}

bool StreamerHepRepTypeTree::addType(HepRepType* type) {
    return true;
}

vector<HepRepType*>* StreamerHepRepTypeTree::getTypes() {
    return NULL;
}

