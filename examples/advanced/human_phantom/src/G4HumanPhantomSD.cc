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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#include "G4HumanPhantomSD.hh"
//#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
//#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4HumanPhantomEnergyDeposit.hh"

G4HumanPhantomSD::G4HumanPhantomSD(G4String name, G4HumanPhantomEnergyDeposit* edepTot)
  :G4VSensitiveDetector(name), energyTotal(edepTot)
{
}

G4HumanPhantomSD::~G4HumanPhantomSD()
{ 
}

void G4HumanPhantomSD::Initialize(G4HCofThisEvent*)
{
}

G4bool G4HumanPhantomSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{  
G4double edep = aStep->GetTotalEnergyDeposit();

  if(edep==0.) return false;
 
  G4String bodypartName = aStep->GetPreStepPoint()->GetTouchable()
   			->GetVolume()->GetLogicalVolume()->GetName();

  energyTotal->Fill(bodypartName, edep);
 return true;
}

void G4HumanPhantomSD::EndOfEvent(G4HCofThisEvent*)
{
}
