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
 * @version $Id: G4HepRepViewer.cc,v 1.7 2002-11-14 05:08:13 duns Exp $
 */

//HepRep
#include "HEPREP/HepRep.h"
#include "HEPREP/HepRepFactory.h"

//G4
#include "G4Types.hh"
//#include "G4VInteractorManager.hh"
#include "G4Scene.hh"
#include "G4HepRep.hh"
#include "G4HepRepSceneHandler.hh"
//This
#include "G4HepRepViewer.hh"

using namespace HEPREP;

G4HepRepViewer::G4HepRepViewer (G4VSceneHandler& scene, const G4String& name)
        : G4VViewer (scene, scene.IncrementViewCount(), name),
          drawn(false) {

#ifdef DEBUG
    G4cout << "G4HepRepViewer::G4HepRepViewer" << G4endl;
#endif
}



G4HepRepViewer::~G4HepRepViewer () {
}


void G4HepRepViewer::ClearView () {
#ifdef DEBUG
    G4cout << "G4HepRepViewer::ClearView" << G4endl;
#endif
}

/***************************************************************************/
// Calculates view representation based on extent of object being
// viewed and (initial) direction of camera.  (Note: it can change
// later due to user interaction via visualization system's GUI.)
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
void G4HepRepViewer::SetView () {
#ifdef DEBUG
    G4cout << "G4HepRepViewer::SetView" << G4endl;
#endif
}



void G4HepRepViewer::DrawView () {
#ifdef DEBUG
    G4cout << "G4HepRepViewer::DrawView" << G4endl;
#endif
    NeedKernelVisit();
    ProcessView();
    drawn = true;
}

void G4HepRepViewer::ShowView () {
#ifdef DEBUG
    G4cout << "G4HepRepViewer::ShowView" << G4endl;
#endif
    G4VViewer::ShowView();

    if (drawn) {
        G4HepRepSceneHandler* sceneHandler = (G4HepRepSceneHandler*)GetSceneHandler();
        sceneHandler->close();
        drawn = false;
    }
}

void G4HepRepViewer::FinishView () {
#ifdef DEBUG
    G4cout << "G4HepRepViewer::FinishView" << G4endl;
#endif
    G4VViewer::FinishView();
}

