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
// $Id: G4HepRepViewer.cc 68043 2013-03-13 14:27:49Z gcosmo $
//

/**
 * @author Mark Donszelmann
 */

//HepRep
#include "HEPREP/HepRep.h"

//G4
#include "G4Types.hh"
//#include "G4VInteractorManager.hh"
#include "G4Scene.hh"
#include "G4HepRep.hh"
#include "G4HepRepMessenger.hh"
#include "G4HepRepSceneHandler.hh"
//This
#include "G4HepRepViewer.hh"

using namespace HEPREP;
using namespace std;

G4HepRepViewer::G4HepRepViewer (G4VSceneHandler& sceneHandler, const G4String& name)
        : G4VViewer (sceneHandler, sceneHandler.IncrementViewCount(), name),
        geometryIncluded(false) {

#ifdef SDEBUG
    cout << "G4HepRepViewer::G4HepRepViewer " << name << endl;
#endif

    // Make changes to view parameters for HepRep...
    fVP.SetCulling(false);
    fDefaultVP.SetCulling(false);
}



G4HepRepViewer::~G4HepRepViewer () {
#ifdef SDEBUG
    cout << "G4HepRepViewer::~G4HepRepViewer" << endl;
#endif
    G4HepRep* pHepRepSystem =
      dynamic_cast<G4HepRep*>(GetSceneHandler()->GetGraphicsSystem());
    if (pHepRepSystem) pHepRepSystem->removeViewer();
}


void G4HepRepViewer::ClearView () {
#ifdef SDEBUG
    cout << "G4HepRepViewer::ClearView" << endl;
#endif
}

void G4HepRepViewer::SetView () {
#ifdef SDEBUG
    cout << "G4HepRepViewer::SetView" << endl;
#endif
}


/* NOTE:
    /run/beamOn         calls ShowView for every event (unless accumulate is set)
    /vis/viewer/flush   calls /vis/viewer/refresh followed by /vis/viewer/update
    /vis/viewer/refresh calls SetView, ClearView, DrawView
    /vis/viewer/update  calls ShowView
*/
void G4HepRepViewer::DrawView () {
#ifdef SDEBUG
    cout << "G4HepRepViewer::DrawView" << endl;
#endif
    if (!geometryIncluded) {
        // draws the geometry
        NeedKernelVisit();
        ProcessView();
        geometryIncluded = true;
    }
}

void G4HepRepViewer::ShowView () {
#ifdef SDEBUG
    cout << "G4HepRepViewer::ShowView" << endl;
#endif
    G4VViewer::ShowView();

    G4HepRepSceneHandler* sceneHandler = dynamic_cast<G4HepRepSceneHandler*>(GetSceneHandler());
    if (sceneHandler) {
      if (sceneHandler->closeHepRep()) {
        sceneHandler->openHepRep();

        G4HepRepMessenger* messenger = G4HepRepMessenger::GetInstance();
        if (messenger->appendGeometry()) geometryIncluded = false;
      }
    }
}

void G4HepRepViewer::FinishView () {
#ifdef SDEBUG
    cout << "G4HepRepViewer::FinishView" << endl;
#endif
    G4VViewer::FinishView();
}

void G4HepRepViewer::reset() {
    geometryIncluded = false;
}  
