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
// G4MSSteppingAction implementation
//
// Author: M.Asai, 5 May 2006
// --------------------------------------------------------------------

#include "G4MSSteppingAction.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Region.hh"
#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"

// --------------------------------------------------------------------
G4MSSteppingAction::G4MSSteppingAction() {}

// --------------------------------------------------------------------
G4MSSteppingAction::~G4MSSteppingAction() {}

// --------------------------------------------------------------------
void G4MSSteppingAction::Initialize(G4bool rSens, G4Region* reg)
{
  regionSensitive = rSens;
  theRegion       = reg;
  length          = 0.;
  x0              = 0.;
  lambda          = 0.;
}

// --------------------------------------------------------------------
void G4MSSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4Region* region =
    preStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetRegion();

  if(regionSensitive && (region != theRegion))
    return;

  G4double stlen       = aStep->GetStepLength();
  const G4Material* material = preStepPoint->GetMaterial();
  length += stlen;
  x0 += stlen / (material->GetRadlen());
  lambda += stlen / (material->GetNuclearInterLength());
}
