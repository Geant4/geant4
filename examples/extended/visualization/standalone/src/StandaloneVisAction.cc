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
// $Id: StandaloneVisAction.cc,v 1.1 2005/10/18 18:09:14 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $

#include "StandaloneVisAction.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Polyhedron.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"

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
