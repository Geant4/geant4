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
/*
 * MoleculeGunMessenger.cc
 *
 *  Created on: 30 janv. 2014
 *      Author: kara
 */

#include "G4MoleculeGunMessenger.hh"
#include "G4MoleculeGun.hh"

#include "G4Tokenizer.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"

//------------------------------------------------------------------------------

G4MoleculeGunMessenger::G4MoleculeGunMessenger(G4MoleculeGun* gun) :
  G4UImessenger("/chem/gun/", "")
{
  fpGunNewGunType = new G4UIcmdWithAString("/chem/gun/newShoot",
                                           this);
  fpMoleculeGun = gun;
}

//------------------------------------------------------------------------------

G4MoleculeGunMessenger::~G4MoleculeGunMessenger()
{
  if (fpGunNewGunType) delete fpGunNewGunType;
}

//------------------------------------------------------------------------------

G4String G4MoleculeGunMessenger::GetCurrentValue(G4UIcommand* /*command*/)
{
  return "";
}

//------------------------------------------------------------------------------

void G4MoleculeGunMessenger::SetNewValue(G4UIcommand* command,
                                         G4String newValue)
{
  if (command == fpGunNewGunType)
  {
    std::istringstream iss (newValue);
    
    G4String shootName;
    iss >> shootName;
    
    G4String shootType;
    iss >> shootType;
        
    if(shootType == "" || shootType.empty())
    {
      CreateNewType<G4Track>(shootName);
    }
    else
    {      
      CreateNewType<G4ContinuousMedium>(shootName);
    }
  }
}

//------------------------------------------------------------------------------

G4MoleculeShootMessenger::G4MoleculeShootMessenger(const G4String& name,
                                                   G4MoleculeGunMessenger*,
                                                   G4shared_ptr<G4MoleculeShoot>
                                                      shoot) :
    G4UImessenger(), fpShoot(shoot)
{
  G4String dir("/chem/gun/");
  dir += name;
  CreateDirectory(dir, "");

  G4String tmp = dir;
  tmp += "/species";
  fpGunSpecies = new G4UIcmdWithAString(tmp, this);

  tmp = dir;
  tmp += "/position";
  fpGunPosition = new G4UIcmdWith3VectorAndUnit(tmp, this);

  tmp = dir;
  tmp += "/time";
  fpGunTime = new G4UIcmdWithADoubleAndUnit(tmp, this);

  tmp = dir;
  tmp += "/number";
  fpGunN = new G4UIcmdWithAnInteger(tmp, this);

  tmp = dir;
  tmp += "/rndmPosition";
  fpGunRdnmPosition = new G4UIcmdWith3VectorAndUnit(tmp, this);

  tmp = dir;
  tmp += "/type";
  fpGunType = new G4UIcmdWithAString(tmp, this);

//  fpShoot.reset(new TG4MoleculeShoot<G4Track>());
}

//------------------------------------------------------------------------------

G4MoleculeShootMessenger::~G4MoleculeShootMessenger()
{
  if (fpGunSpecies) delete fpGunSpecies;
  if (fpGunPosition) delete fpGunPosition;
  if (fpGunTime) delete fpGunTime;
  if (fpGunN) delete fpGunN;
}

//------------------------------------------------------------------------------

void G4MoleculeShootMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fpGunSpecies)
  {
    fpShoot->fMoleculeName = newValue;
  }
  else if (command == fpGunPosition)
  {
    fpShoot->fPosition = fpGunPosition->GetNew3VectorValue(newValue);
  }
  else if(command == fpGunRdnmPosition)
  {
    fpShoot->fBoxSize = new G4ThreeVector(fpGunRdnmPosition->GetNew3VectorValue(newValue));
  }
  else if (command == fpGunTime)
  {
    fpShoot->fTime = fpGunTime->GetNewDoubleValue(newValue);
  }
  else if (command == fpGunN)
  {
    fpShoot->fNumber = fpGunN->GetNewIntValue(newValue);
  }
  else if (command == fpGunType)
  {
    if(newValue == "CM")
    {
//      G4cout << "**** Change type" << G4endl;
//      TG4MoleculeShoot<G4ContinuousMedium>* casted = static_cast<TG4MoleculeShoot<G4ContinuousMedium>*>(fpShoot.get());
//      fpShoot.reset(casted);
      fpShoot = fpShoot.get()->ChangeType<G4ContinuousMedium>();
    }
  }
}

//------------------------------------------------------------------------------

G4String G4MoleculeShootMessenger::GetCurrentValue(G4UIcommand* command)
{
  if (command == fpGunSpecies)
  {
    return fpShoot->fMoleculeName;
  }
  else if (command == fpGunPosition)
  {
    return fpGunPosition->ConvertToStringWithBestUnit(fpShoot->fPosition);
  }
  else if (command == fpGunRdnmPosition)
  {
    if(fpShoot->fBoxSize)
    {
      return fpGunRdnmPosition->ConvertToStringWithBestUnit(*fpShoot->fBoxSize);
    }
    return fpGunRdnmPosition->ConvertToStringWithBestUnit(G4ThreeVector());
  }
  else if (command == fpGunTime)
  {
    return fpGunTime->ConvertToStringWithBestUnit(fpShoot->fTime);
  }
  else if (command == fpGunN)
  {
    return fpGunN->ConvertToString(fpShoot->fNumber);
  }
  return "";
}

//------------------------------------------------------------------------------

