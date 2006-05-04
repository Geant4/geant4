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
// $Id: G4MSSteppingAction.cc,v 1.1 2006-05-04 19:42:47 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//


#include "G4MSSteppingAction.hh"

#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Region.hh"
#include "G4Material.hh"

G4MSSteppingAction::G4MSSteppingAction()
{
  Initialize(false,0);
}

G4MSSteppingAction::~G4MSSteppingAction()
{;}
 
void G4MSSteppingAction::Initialize(G4bool rSens,G4Region* reg)
{
  regionSensitive = rSens;
  theRegion = reg;
  length = 0.;
  x0 = 0.;
  lambda = 0.;
}
  
void G4MSSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4Region* region = preStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetRegion();

  if(regionSensitive && (region!=theRegion)) return;

  G4double stlen = aStep->GetStepLength();
  G4Material* material = preStepPoint->GetMaterial();
  length += stlen;
  x0 += stlen/(material->GetRadlen());
  lambda += stlen/(material->GetNuclearInterLength());
}

