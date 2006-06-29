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
// $Id: Tst01ParticleGunMessenger.cc,v 1.3 2006-06-29 19:32:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst01ParticleGunMessenger.hh"
#include "Tst01ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "G4ios.hh"


Tst01ParticleGunMessenger::Tst01ParticleGunMessenger(Tst01ParticleGun * fPtclGun)
  :fParticleGun(fPtclGun)
{
  particleTable = G4ParticleTable::GetParticleTable();

  properTimeCmd= new G4UIcmdWithADoubleAndUnit("/gun/decayTime",this);
  properTimeCmd->SetGuidance("Set pre-assigned decay time.");
  properTimeCmd->SetParameterName("t",true,true);
  properTimeCmd->SetDefaultUnit("ns");
  properTimeCmd->SetRange("t>=0.");
}

Tst01ParticleGunMessenger::~Tst01ParticleGunMessenger()
{
  delete properTimeCmd;
}

void Tst01ParticleGunMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command==properTimeCmd ){
    fParticleGun->SetDecayProperTime(properTimeCmd->GetNewDoubleValue(newValues)); 
  }
}

G4String Tst01ParticleGunMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==properTimeCmd ){
    cv = properTimeCmd ->ConvertToString(fParticleGun->GetDecayProperTime(),"ns"); 
  }    
  return cv;
}



