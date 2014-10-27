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

G4MoleculeGunMessenger::G4MoleculeGunMessenger()
{
  fpGunDir = new G4UIdirectory("/process/em/dna/chem/gun/");
  fpGunNewGunType = new G4UIcmdWithAString("/process/em/dna/chem/gun/newType",
                                           this);
}

G4MoleculeGunMessenger::~G4MoleculeGunMessenger()
{
  if (fpGunDir) delete fpGunDir;
  if (fpGunNewGunType) delete fpGunNewGunType;
}

G4String G4MoleculeGunMessenger::GetCurrentValue(G4UIcommand* /*command*/)
{
  return "";
}

G4MoleculeGunMessenger::MultipleGun* G4MoleculeGunMessenger::CreateNewType(const G4String& name)
{
  MultipleGun* multiGun = new MultipleGun(name, this);
  fMultipleGun.push_back(multiGun);
  return multiGun;
}

void G4MoleculeGunMessenger::SetNewValue(G4UIcommand* command,
                                         G4String newValue)
{
  if (command == fpGunNewGunType)
  {
    CreateNewType(newValue);
  }
}

G4MoleculeGunMessenger::MultipleGun::MultipleGun(const G4String& name,
                                                 G4MoleculeGunMessenger*)
{
  G4String dir("/process/em/dna/chem/gun/");
  dir += name;
  fpGunType = new G4UIdirectory(name);

  G4String tmp = dir;
  tmp += "/moleculeModel";
  fpGunMoleculeModel = new G4UIcmdWithAString(tmp, this);
  tmp = dir;
  tmp += "/position";
  fpGunPosition = new G4UIcmdWith3VectorAndUnit(tmp, this);
  tmp = dir;
  tmp += "/time";
  fpGunTime = new G4UIcmdWithADoubleAndUnit(tmp, this);
  tmp = dir;
  tmp += "/number";
  fpGunN = new G4UIcmdWithAnInteger(tmp, this);

  fMoleculeName = "";
  fTime = 0;
  fNumber = 0;
}

G4MoleculeGunMessenger::MultipleGun::~MultipleGun()
{
  if (fpGunMoleculeModel) delete fpGunMoleculeModel;
  if (fpGunPosition) delete fpGunPosition;
  if (fpGunTime) delete fpGunTime;
  if (fpGunN) delete fpGunN;
}

void G4MoleculeGunMessenger::MultipleGun::SetNewValue(G4UIcommand* command,
                                                      G4String newValue)
{

  if (command == fpGunMoleculeModel)
  {
    fMoleculeName = newValue;
  }
  else if (command == fpGunPosition)
  {
    fPosition = fpGunPosition->GetNew3VectorValue(newValue);
  }
  else if (command == fpGunTime)
  {
    fTime = fpGunTime->GetNewDoubleValue(newValue);
  }
  else if (command == fpGunN)
  {
    fNumber = fpGunN->GetNewIntValue(newValue);
  }
}

G4String G4MoleculeGunMessenger::MultipleGun::GetCurrentValue(G4UIcommand* command)
{

  if (command == fpGunMoleculeModel)
  {
    return fMoleculeName;
  }
  else if (command == fpGunPosition)
  {
    return fpGunPosition->ConvertToStringWithBestUnit(fPosition);
  }
  else if (command == fpGunTime)
  {
    return fpGunTime->ConvertToStringWithBestUnit(fTime);
  }
  else if (command == fpGunN)
  {
    return fpGunN->ConvertToString(fNumber);
  }
  return "";
}

void G4MoleculeGunMessenger::DefineTracks(G4MoleculeGun* gun)
{
  for (size_t i = 0; i < fMultipleGun.size(); i++)
  {
    fMultipleGun[i]->DefineTracks(gun);
  }
}

void G4MoleculeGunMessenger::MultipleGun::DefineTracks(G4MoleculeGun* gun)
{
  gun->AddNMolecules(fNumber, fMoleculeName, fPosition, fTime);
}
