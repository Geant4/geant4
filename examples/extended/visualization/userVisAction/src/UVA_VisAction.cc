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
// $Id: UVA_VisAction.cc,v 1.3 2006-06-29 17:47:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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
    G4VisAttributes va1(G4Colour::Red());
    va1.SetForceSolid(true);
    pVisManager->Draw(logo,va1,G4Translate3D(0,-1*m,4.5*m));

    G4Text text("My beautiful logo");
    G4VisAttributes va2(G4Colour::Magenta());
    text.SetVisAttributes(va2);
    text.SetScreenSize(12.);
    pVisManager->Draw(text,G4Translate3D(0,0,3.5*m));

  }
}
