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

#include "StreamerHepRep.h"

using namespace std;
using namespace HEPREP;

StreamerHepRep::StreamerHepRep(HepRepWriter* streamer) {
    streamer->write(this);
}

StreamerHepRep::~StreamerHepRep() {
}

HepRep* StreamerHepRep::copy(HepRepSelectFilter*) {
    return NULL;
}

vector<string>* StreamerHepRep::getLayerOrder() {
    return &layers;
}

bool StreamerHepRep::addLayer(string layer) {
    layers.push_back(layer);
    return true;
}

bool StreamerHepRep::addTypeTree(HepRepTypeTree*) {
    return true;
}

void StreamerHepRep::removeTypeTree(HepRepTypeTree*) {
}

HepRepTypeTree* StreamerHepRep::getTypeTree(string, string) {
    return NULL;
}

vector<HepRepTypeTree*>* StreamerHepRep::getTypeTrees() {
    return NULL;
}

HepRepType* StreamerHepRep::getType(string) {
    return NULL;
}

bool StreamerHepRep::addInstanceTree(HepRepInstanceTree*) {
    return true;
}

void StreamerHepRep::removeInstanceTree(HepRepInstanceTree*) {
}

HepRepInstanceTree* StreamerHepRep::getInstanceTreeTop(string, string) {
    return NULL;
}

HepRepInstanceTree* StreamerHepRep::getInstances(string, string,
                                       vector<string>) {
    return NULL;
}

HepRepInstanceTree* StreamerHepRep::getInstancesAfterAction(
                                string,
                                string,
                                vector<string>,
                                vector<HepRepAction*>,
                                bool,
                                bool,
                                bool,
                                vector<string>) {
    return NULL;
}

string StreamerHepRep::checkForException() {
    return NULL;
}

vector<HepRepInstanceTree*>* StreamerHepRep::getInstanceTrees() {
    return NULL;
}

