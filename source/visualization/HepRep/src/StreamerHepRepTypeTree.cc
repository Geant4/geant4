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
//

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

