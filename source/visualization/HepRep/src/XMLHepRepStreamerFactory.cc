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

#include "XMLHepRepStreamerFactory.h"
#include "XMLHepRepStreamer.h"

#include "StreamerHepRepPoint.h"
#include "StreamerHepRepInstance.h"
#include "StreamerHepRepInstanceTree.h"
#include "StreamerHepRepType.h"
#include "StreamerHepRepTypeTree.h"
#include "StreamerHepRep.h"
#include "DefaultHepRepAction.h"
#include "DefaultHepRepTreeID.h"

using namespace std;
using namespace HEPREP;


XMLHepRepStreamerFactory::XMLHepRepStreamerFactory() {
}

XMLHepRepStreamerFactory::~XMLHepRepStreamerFactory() {
}

HepRepWriter* XMLHepRepStreamerFactory::createHepRepWriter(ostream* out) {
    this->streamer = new XMLHepRepStreamer(out);
    return streamer;
}

HepRepPoint* XMLHepRepStreamerFactory::createHepRepPoint (HepRepInstance* instance,
                               double x, double y, double z) {
    return new StreamerHepRepPoint(streamer, instance, x, y, z);
}

HepRepInstance* XMLHepRepStreamerFactory::createHepRepInstance (HepRepInstance* parent, HepRepType* type) {
    return new StreamerHepRepInstance(streamer, parent, type);
}

HepRepInstance* XMLHepRepStreamerFactory::createHepRepInstance (HepRepInstanceTree* parent, HepRepType* type) {
    return new StreamerHepRepInstance(streamer, parent, type);
}

HepRepTreeID* XMLHepRepStreamerFactory::createHepRepTreeID (string name, string version, string qualifier) {
    return new DefaultHepRepTreeID(name, version, qualifier);
}

HepRepAction* XMLHepRepStreamerFactory::createHepRepAction (string name, string expression) {
    return new DefaultHepRepAction(name, expression);
}

HepRepInstanceTree* XMLHepRepStreamerFactory::createHepRepInstanceTree (string name, string version,
                                                    HepRepTreeID* typeTreeID) {
    return new StreamerHepRepInstanceTree(streamer, name, version, typeTreeID);
}

HepRepType* XMLHepRepStreamerFactory::createHepRepType (HepRepType* parent, string name) {
    return new StreamerHepRepType(streamer, parent, name);
}

HepRepTypeTree* XMLHepRepStreamerFactory::createHepRepTypeTree (HepRepTreeID* treeID) {
    return new StreamerHepRepTypeTree(streamer, treeID);
}

HepRep* XMLHepRepStreamerFactory::createHepRep () {
    return new StreamerHepRep(streamer);
}

