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
// $Id: G4VUserDetectorConstruction.cc,v 1.4 2006-04-26 15:24:24 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserPhysicalVolume.hh"

G4VUserDetectorConstruction::G4VUserDetectorConstruction()
{;}

G4VUserDetectorConstruction::~G4VUserDetectorConstruction()
{;}

void G4VUserDetectorConstruction::RegisterParallelWorld(G4VUserParallelWorld* aPW)
{
  std::vector<G4VUserParallelWorld*>::iterator pwItr;
  for(pwIte=parallelWorld.begin();pwItr!=parallelWorld.end();pwItr++)
  {
    if((*pwItr)->GetName()==aPW->GetName())
    {
      G4String eM = "A parallel world <";
      eM += aPW->GetName();
      eM += "> is already registered to the user detector construction.";
      G4Exception("G4VUserDetectorConstruction::RegisterParallelWorld",
                  "RunUDet000",FatalErrorInArgument,eM);
    }
  }
  parallelWorld.push_back(aPW);
}

void G4VUserDetectorConstruction::ConstructParallelGeometries()
{
  std::vector<G4VUserParallelWorld*>::iterator pwItr;
  for(pwIte=parallelWorld.begin();pwItr!=parallelWorld.end();pwItr++)
  {
    pwItr->Construct();
  }
}

