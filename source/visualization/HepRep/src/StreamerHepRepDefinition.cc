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

StreamerHepRepDefinition::StreamerHepRepDefinition(HepRepWriter* streamer)
    : StreamerHepRepAttribute(streamer), streamer(streamer) {
}

StreamerHepRepDefinition::~StreamerHepRepDefinition() {
}

vector<HepRepAttDef *>* StreamerHepRepDefinition::getAttDefsFromNode() {
    return NULL;
}

bool StreamerHepRepDefinition::addAttDef(HepRepAttDef* hepRepAttDef) {
    streamer->write(hepRepAttDef);
    delete hepRepAttDef;
    return true;
}

bool StreamerHepRepDefinition::addAttDef(string name, string desc, string type, string extra) {
    addAttDef(new DefaultHepRepAttDef(name, desc, type, extra));
    return true;
}

HepRepAttDef* StreamerHepRepDefinition::getAttDefFromNode(string lowerCaseName) {
    return NULL;
}

HepRepAttDef* StreamerHepRepDefinition::getAttDef(string name) {
    return NULL;
}

