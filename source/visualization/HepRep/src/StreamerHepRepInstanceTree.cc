
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

HepRepInstanceTree* StreamerHepRepInstanceTree::copy(HepRep* heprep, HepRepSelectFilter* filter) {
    return NULL;
}

bool StreamerHepRepInstanceTree::addInstance(HepRepInstance* instance) {
    return true;
}

void StreamerHepRepInstanceTree::removeInstance(HepRepInstance* instance) {
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

