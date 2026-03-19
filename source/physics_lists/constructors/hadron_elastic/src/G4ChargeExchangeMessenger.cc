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
//---------------------------------------------------------------------------
//
// ClassName:   G4ChargeExchangeMessenger
//
// Author: 2023 V. Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4ChargeExchangeMessenger.hh"
#include "G4ChargeExchangePhysics.hh"

G4ChargeExchangeMessenger::G4ChargeExchangeMessenger(G4ChargeExchangePhysics* a)
  : theB(a)
{
  // general stuff.
  aDir = new G4UIdirectory("/physics_lists/cex/", false);
  aDir->SetGuidance("tailoring the hadronic charge exchange processes.");

  fCmd = new G4UIcmdWithADouble("/physics_lists/cex/PionBiasFactor",this);
  fCmd->SetGuidance("Pion charge exchange cross section factor");
  fCmd->AvailableForStates(G4State_PreInit);
  fCmd->SetToBeBroadcasted(false);

  fCmd1 = new G4UIcmdWithADouble("/physics_lists/cex/KaonBiasFactor",this);
  fCmd1->SetGuidance("Kaon charge exchange cross section factor");
  fCmd1->AvailableForStates(G4State_PreInit);
  fCmd1->SetToBeBroadcasted(false);
  
  lCmd = new G4UIcmdWithADoubleAndUnit("/process/cex/LowEnergyLimit",this);
  lCmd->SetGuidance("Low-energy energy limit for charge exchange process");
  lCmd->SetParameterName("cexLowE",true);
  lCmd->SetUnitCategory("Energy");
  lCmd->AvailableForStates(G4State_PreInit);
  lCmd->SetToBeBroadcasted(false);
}

G4ChargeExchangeMessenger::~G4ChargeExchangeMessenger()
{
  delete fCmd;
  delete fCmd1;
  delete lCmd;
  delete aDir;
}

void G4ChargeExchangeMessenger::SetNewValue(G4UIcommand* aComm, G4String aS)
{
  if (aComm == fCmd) {
    theB->SetPionCrossSectionFactor(fCmd->GetNewDoubleValue(aS));
  } else if (aComm == fCmd1) {
    theB->SetKaonCrossSectionFactor(fCmd1->GetNewDoubleValue(aS));
  } else if (aComm == lCmd) {
    theB->SetLowEnergyLimit(lCmd->GetNewDoubleValue(aS));
  }
}
