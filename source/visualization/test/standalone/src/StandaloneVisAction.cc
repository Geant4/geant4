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
// $Id: StandaloneVisAction.cc,v 1.1 2005-02-19 21:50:20 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "StandaloneVisAction.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Polyhedron.hh"
#include "G4Box.hh"

void StandaloneVisAction::Draw() {
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if (pVisManager) {
    pVisManager->Draw(G4Box("box",2*m,2*m,2*m),
		      G4VisAttributes(G4Colour(1,1,0)));
    G4Polyhedron* polyhedron = G4Box("box",2*m,2*m,2*m).CreatePolyhedron();
    polyhedron = new G4Polyhedron(polyhedron->subtract(*polyhedron));
  }
}
