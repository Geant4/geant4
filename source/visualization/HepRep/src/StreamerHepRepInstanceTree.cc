
#include "StreamerHepRepInstanceTree.h"

using namespace std;
using namespace HEPREP;

StreamerHepRepInstanceTree::StreamerHepRepInstanceTree(HepRepWriter* stream, string name, string version, HepRepTreeID* typeTree)
    : DefaultHepRepTreeID(name, version), streamer(stream), typeTree(typeTree) {

    stream->write(this);
}

StreamerHepRepInstanceTree::~StreamerHepRepInstanceTree() {
}

HepRepTreeID* StreamerHepRepInstanceTree::copy() {
    return DefaultHepRepTreeID::copy();
}

HepRepInstanceTree* StreamerHepRepInstanceTree::copy(HepRep*, HepRepSelectFilter*) {
    return NULL;
}

bool StreamerHepRepInstanceTree::addInstance(HepRepInstance*) {
    return true;
}

void StreamerHepRepInstanceTree::removeInstance(HepRepInstance*) {
}

vector<HepRepInstance*>* StreamerHepRepInstanceTree::getInstances() {
    return NULL;
}

bool StreamerHepRepInstanceTree::addInstanceTree(HepRepTreeID* treeID) {
    streamer->write(treeID);
    return true;
}

HepRepTreeID* StreamerHepRepInstanceTree::getTypeTree() {
    return typeTree;
}

vector<HepRepInstanceTree*>* StreamerHepRepInstanceTree::getInstanceTrees() {
    return NULL;
}

