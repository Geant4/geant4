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
//

#include "G4RayTracerQtViewer.hh"

#include "G4VSceneHandler.hh"
#include "G4Scene.hh"
#include "G4TheMTRayTracer.hh"
#include "G4RTJpegMaker.hh"
#include "G4RTQtScanner.hh"
#include "G4UImanager.hh"

#include <cstdlib>

G4RayTracerQtViewer::G4RayTracerQtViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
G4RayTracerViewer
(sceneHandler, name,
 G4TheMTRayTracer::Instance(new G4RTJpegMaker, new G4RTQtScanner))
{}

G4RayTracerQtViewer::~G4RayTracerQtViewer() {}

void G4RayTracerQtViewer::Initialise() {

  G4RayTracerViewer::Initialise();

  fVP.SetAutoRefresh(true);
  fDefaultVP.SetAutoRefresh(true);

  // Set up Qt Window...
  G4RTQtScanner* theQtScanner = (G4RTQtScanner*)theTracer->GetScanner();
  if (!theQtScanner->GetQtWindow(fName,fVP)) {
    G4cerr << "G4RayTracerQtViewer::Initialise: No scanner" << G4endl;
    fViewId = -1;  // This flags an error.
    return;
  }
}
