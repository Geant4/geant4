
#include "StreamerHepRepType.h"

using namespace std;
using namespace HEPREP;

StreamerHepRepType::StreamerHepRepType(HepRepWriter* stream, HepRepType* parent, string name)
    : StreamerHepRepDefinition(stream), parent(parent), name(name) {
    this->description = "No Description";
    this->infoURL = "No Info URL";

    stream->write(this);
}

StreamerHepRepType::~StreamerHepRepType() {
}

HepRepType* StreamerHepRepType::getSuperType() {
    return parent;
}

HepRepAttDef* StreamerHepRepType::getAttDef(string key) {
    return NULL;
}

/**
 * searched for a value with given name. Search up the type tree if needed.
 */
HepRepAttValue* StreamerHepRepType::getAttValue(string key) {
    return NULL;
}

HepRepType* StreamerHepRepType::copy(HepRep* heprep, HepRepType* type) {
    return NULL;
}

string StreamerHepRepType::getName() {
    return name;
}

string StreamerHepRepType::getDescription() {
    return description;
}

void StreamerHepRepType::setDescription(string desc) {
    this->description = desc;
}

string StreamerHepRepType::getInfoURL() {
    return infoURL;
}

void StreamerHepRepType::setInfoURL(string info) {
    this->infoURL = info;
}

bool StreamerHepRepType::addType(HepRepType* type) {
    return true;
}

vector<HepRepType*>* StreamerHepRepType::getTypes() {
    return NULL;
}

