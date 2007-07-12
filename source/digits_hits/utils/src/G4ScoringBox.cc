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
// $Id: G4ScoringBox.cc,v 1.1 2007-07-12 05:16:46 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4ScoringBox.hh"
#include "G4VPhysicalVolume.hh"

G4ScoringBox::G4ScoringBox(G4String wName)
:G4VScoringMesh(wName)
{
  fShape = boxMesh;
}

G4ScoringBox::~G4ScoringBox()
{
}

void G4ScoringBox::Construct(G4VPhysicalVolume* fWorldPhys)
{
  if(fConstructed) 
  {
    G4cerr << fWorldPhys->GetName() << G4endl;
    G4Exception(fWorldName+" has already been built.");
  }
  fConstructed = true;

}



