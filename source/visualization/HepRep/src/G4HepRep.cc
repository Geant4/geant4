//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4HepRep.cc 78838 2014-01-28 08:46:17Z gcosmo $
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
                             G4VGraphicsSystem::fileWriter),
          sceneHandler(NULL),
          viewer(NULL) {
		G4HepRepMessenger::GetInstance();
}

G4HepRep::~G4HepRep () {
}

G4VSceneHandler* G4HepRep::CreateSceneHandler (const G4String& name) {
    if (sceneHandler != NULL) {
        cout << "G4HepRep::CreateSceneHandler: Cannot create more than one G4HepRepSceneHandler" << endl;
        return NULL;
    }
    sceneHandler = new G4HepRepSceneHandler (*this, name);
    return sceneHandler;
}

G4VViewer* G4HepRep::CreateViewer (G4VSceneHandler& scene, const G4String& name) {
    if (viewer != NULL) {
        cout << "G4HepRep::CreateViewer: Cannot create more than one G4HepRepViewer" << endl;
        return NULL;
    }
    viewer  = new G4HepRepViewer ((G4HepRepSceneHandler&)scene, name);
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

