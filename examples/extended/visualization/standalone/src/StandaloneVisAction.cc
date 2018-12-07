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
/// \file visualization/standalone/src/StandaloneVisAction.cc
/// \brief Implementation of the StandaloneVisAction class
//
//

#include "StandaloneVisAction.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Polyhedron.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StandaloneVisAction::Draw() {
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if (pVisManager) {

    // Simple box...
    pVisManager->Draw(G4Box("box",2*m,2*m,2*m),
                      G4VisAttributes(G4Colour(1,1,0)));

    // Boolean solid...
    G4Box boxA("boxA",3*m,3*m,3*m);
    G4Box boxB("boxB",1*m,1*m,1*m);
    G4SubtractionSolid subtracted("subtracted_boxes",&boxA,&boxB,
                       G4Translate3D(3*m,3*m,3*m));
    pVisManager->Draw(subtracted,
                      G4VisAttributes(G4Colour(0,1,1)),
                      G4Translate3D(-6*m,-6*m,-6*m));

    // Same, but explicit polyhedron...
    G4Polyhedron* pA = G4Box("boxA",3*m,3*m,3*m).CreatePolyhedron();
    G4Polyhedron* pB = G4Box("boxB",1*m,1*m,1*m).CreatePolyhedron();
    pB->Transform(G4Translate3D(3*m,3*m,3*m));
    G4Polyhedron* pSubtracted = new G4Polyhedron(pA->subtract(*pB));
    G4VisAttributes subVisAtts(G4Colour(0,1,1));
    pSubtracted->SetVisAttributes(&subVisAtts);
    pVisManager->Draw(*pSubtracted,G4Translate3D(6*m,6*m,6*m));
    delete pA;
    delete pB;
    delete pSubtracted;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
