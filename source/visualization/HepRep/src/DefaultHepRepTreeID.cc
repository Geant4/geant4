//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
#include "DefaultHepRepTreeID.h"

using namespace std;
using namespace HEPREP;

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

HepRepTreeID* DefaultHepRepTreeID::copy() {
    return NULL;
}

