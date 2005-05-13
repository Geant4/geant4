
#include <iostream>

#include "cheprep/DefaultHepRepTreeID.h"

using namespace std;
using namespace HEPREP;

namespace cheprep {

DefaultHepRepTreeID::DefaultHepRepTreeID(string name, string version, string qualifier)
    : name(name), version(version), qualifier(qualifier) {
}

DefaultHepRepTreeID::~DefaultHepRepTreeID() {
}

string DefaultHepRepTreeID::getQualifier() {
    return qualifier;
}

void DefaultHepRepTreeID::setQualifier(string qual) {
    this->qualifier = qual;
}

string DefaultHepRepTreeID::getName() {
    return name;
}

string DefaultHepRepTreeID::getVersion() {
    return version;
}

} // cheprep
