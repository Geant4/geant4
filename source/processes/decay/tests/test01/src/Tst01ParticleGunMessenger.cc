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
// $Id: Tst01ParticleGunMessenger.cc,v 1.2 2001-07-11 10:02:28 gunter Exp $
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



