// Copyright FreeHEP, 2005.

#include <iostream>

#include "cheprep/DefaultHepRepTreeID.h"

using namespace std;
using namespace HEPREP;

/**
 * @author Mark Donszelmann
 */
namespace cheprep {

DefaultHepRepTreeID::DefaultHepRepTreeID(string aName, string aVersion, string aQualifier)
    : name(aName), version(aVersion), qualifier(aQualifier) {
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
