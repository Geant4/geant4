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
//

/**
 * @author Mark Donszelmann
 */

//G4
#include "G4Types.hh"
#include "G4HepRepViewer.hh"
#include "G4HepRepSceneHandler.hh"

//This
#include "G4HepRep.hh"

using namespace HEPREP;

G4HepRep::G4HepRep ()
        : G4VGraphicsSystem ("G4HepRep",
			     "HepRepXML",
        "HepRep Generic Driver for XML, RMI and CORBA",
        G4VGraphicsSystem::threeD) {
}

G4HepRep::~G4HepRep () {
}

G4VSceneHandler* G4HepRep::CreateSceneHandler (const G4String& name) {
    G4HepRepSceneHandler* scene = new G4HepRepSceneHandler (*this, name);
    return scene;
}

G4VViewer* G4HepRep::CreateViewer (G4VSceneHandler& scene, const G4String& name) {
    G4VViewer* view  = new G4HepRepViewer ((G4HepRepSceneHandler&)scene, name);
    if (view) {
        if (view->GetViewId() < 0) {
            G4cout << "G4HepRep::CreateViewer: ERROR flagged by negative view ID." << G4endl;
            delete view;
            view = NULL;
        }
    } else {
        G4cout << "G4HepRep::CreateViewer: null pointer on new G4HepRepViewer." << G4endl;
    }
    return view;
}
