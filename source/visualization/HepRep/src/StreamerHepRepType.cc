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

