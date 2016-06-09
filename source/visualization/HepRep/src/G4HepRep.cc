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
// $Id: G4HepRep.cc,v 1.21 2004/05/27 05:55:20 duns Exp $
// GEANT4 tag $Name: geant4-06-02 $
//

/**
 * We allow only one HepRep Scene handler and one attached Viewer.
 * The sceneHandler associates with the output file
 *
 * @author Mark Donszelmann
 */

//G4
#include "G4Types.hh"
#include "G4HepRepViewer.hh"
#include "G4HepRepSceneHandler.hh"
#include "G4HepRepMessenger.hh"

//This
#include "G4HepRep.hh"

using namespace HEPREP;
using namespace std;

G4HepRep::G4HepRep ()
        : G4VGraphicsSystem ("G4HepRep",
                             "HepRepXML",
        	             "HepRep Generic Driver for XML, RMI and CORBA",
                             G4VGraphicsSystem::threeD),
          sceneHandler(NULL),
          viewer(NULL) {
    messenger = new G4HepRepMessenger();
}

G4HepRep::~G4HepRep () {
    delete messenger;
}

G4VSceneHandler* G4HepRep::CreateSceneHandler (const G4String& name) {
    if (sceneHandler != NULL) {
        cout << "G4HepRep::CreateSceneHandler: Cannot create more than one G4HepRepSceneHandler" << endl;
        return NULL;
    }
    sceneHandler = new G4HepRepSceneHandler (*this, *messenger, name);
    return sceneHandler;
}

G4VViewer* G4HepRep::CreateViewer (G4VSceneHandler& scene, const G4String& name) {
    if (viewer != NULL) {
        cout << "G4HepRep::CreateViewer: Cannot create more than one G4HepRepViewer" << endl;
        return NULL;
    }
    viewer  = new G4HepRepViewer ((G4HepRepSceneHandler&)scene, *messenger, name);
    return viewer;
}

void G4HepRep::removeSceneHandler() {
    // actual deletion is done in VisManager
    sceneHandler = NULL;
}

void G4HepRep::removeViewer() {
    // actual deletion is done in VisManager
    viewer = NULL;
}

