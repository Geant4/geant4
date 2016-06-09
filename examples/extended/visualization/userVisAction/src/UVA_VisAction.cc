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
// $Id: UVA_VisAction.cc,v 1.1 2005/10/18 18:09:14 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $

#include "UVA_VisAction.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4Text.hh"

void UVA_VisAction::Draw() {
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if (pVisManager) {

    // A simple logo...
    G4Orb orb("my_logo_orb", 1*m);
    G4Box box("my_cut_box", 1*m, 1*m, 1*m);
    G4SubtractionSolid logo("my_logo", &orb, &box, G4Translate3D(-1*m,1*m,1*m));
    G4VisAttributes va1(G4Colour::Red);
    va1.SetForceSolid(true);
    pVisManager->Draw(logo,va1,G4Translate3D(0,-1*m,4.5*m));

    G4Text text("My beautiful logo");
    G4VisAttributes va2(G4Colour::Magenta);
    text.SetVisAttributes(va2);
    text.SetScreenSize(12.);
    pVisManager->Draw(text,G4Translate3D(0,0,3.5*m));

  }
}
