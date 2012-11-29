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
// $Id$
//

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserParallelWorld.hh"

G4VUserDetectorConstruction::G4VUserDetectorConstruction()
{;}

G4VUserDetectorConstruction::~G4VUserDetectorConstruction()
{;}

void G4VUserDetectorConstruction::RegisterParallelWorld(G4VUserParallelWorld* aPW)
{
  std::vector<G4VUserParallelWorld*>::iterator pwItr;
  for(pwItr=parallelWorld.begin();pwItr!=parallelWorld.end();pwItr++)
  {
    if((*pwItr)->GetName()==aPW->GetName())
    {
      G4String eM = "A parallel world <";
      eM += aPW->GetName();
      eM += "> is already registered to the user detector construction.";
      G4Exception("G4VUserDetectorConstruction::RegisterParallelWorld",
                  "Run0051",FatalErrorInArgument,eM);
    }
  }
  parallelWorld.push_back(aPW);
}

G4int G4VUserDetectorConstruction::ConstructParallelGeometries()
{
  G4int nP = 0;
  std::vector<G4VUserParallelWorld*>::iterator pwItr;
  for(pwItr=parallelWorld.begin();pwItr!=parallelWorld.end();pwItr++)
  {
    (*pwItr)->Construct();
    nP++;
  }
  return nP;
}

G4int G4VUserDetectorConstruction::GetNumberOfParallelWorld() const
{ return parallelWorld.size(); }

G4VUserParallelWorld* G4VUserDetectorConstruction::GetParallelWorld(G4int i) const
{
  if(i<0||i>=GetNumberOfParallelWorld()) return 0;
  return parallelWorld[i];
}

