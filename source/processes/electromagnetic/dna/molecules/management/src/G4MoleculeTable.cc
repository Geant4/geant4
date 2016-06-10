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
 * G4MoleculeTable.cc
 *
 *  Created on: 23 oct. 2013
 *      Author: kara
 */

#include "G4MoleculeTable.hh"

G4MoleculeTable* G4MoleculeTable::fpgMoleculeTable(0);

G4MoleculeTable::G4MoleculeTable()
{
  // TODO Auto-generated constructor stub

}

G4MoleculeTable::~G4MoleculeTable()
{
  // TODO Auto-generated destructor stub
}

G4MoleculeTable* G4MoleculeTable::Instance()
{
  if (!fpgMoleculeTable) fpgMoleculeTable = new G4MoleculeTable;
  return fpgMoleculeTable;
}

G4MoleculeTable* G4MoleculeTable::GetMoleculeTable()
{
  return Instance();
}

G4MoleculeDefinition* G4MoleculeTable::CreateMoleculeDefinition(const G4String& name,
                                                                double diffusion_coefficient)
{
  MoleculeDefTable::iterator it = fMoleculeDefTable.find(name);
  G4MoleculeDefinition* definition(0);
  if (it == fMoleculeDefTable.end())
  {
    definition = new G4MoleculeDefinition(name, -1 /* mass*/,
                                          diffusion_coefficient);
    fMoleculeDefTable[name] = definition;
  }
  else
  {
    // exception
    G4ExceptionDescription description;
    description << "The molecule definition " << name
                << " was already recorded in the table" << G4endl;
    G4Exception("G4MoleculeTable::CreateMoleculeDefinition",
                "DEFINITION_ALREADY_CREATED", FatalException, description);
  }
  return definition;
}

G4Molecule* G4MoleculeTable::CreateMoleculeModel(const G4String& name,
                                                 G4MoleculeDefinition* molDef,
                                                 int charge,
                                                 double diffusion_coefficient)
{

  G4Molecule* molecule = new G4Molecule(molDef, charge);
  if (diffusion_coefficient != -1) molecule->SetDiffusionCoefficient(
      diffusion_coefficient);

  RecordMoleculeModel(name, molecule);

  return molecule;
}

G4Molecule* G4MoleculeTable::CreateMoleculeModel(const G4String& name,
                                                 G4MoleculeDefinition* molDef)
{

  G4Molecule* molecule = new G4Molecule(molDef);
  RecordMoleculeModel(name, molecule);

  return molecule;
}

G4MoleculeDefinition* G4MoleculeTable::GetMoleculeDefinition(const G4String& name, bool mustExist)
{
  MoleculeDefTable::iterator it = fMoleculeDefTable.find(name);
  G4MoleculeDefinition* definition(0);
  if (it != fMoleculeDefTable.end())
  {
    definition = it->second;
  }
  else if(mustExist)
  {
    // exception
    G4ExceptionDescription description;
    description << "The molecule defintion " << name
                << " was NOT recorded in the table" << G4endl;
    G4Exception("G4MoleculeTable::CreateMoleculeModel", "MOLECULE_DEFINITION_NOT_CREATED",
                FatalException, description);
  }
  return definition;
}

void G4MoleculeTable::RecordMoleculeModel(const G4String& name,
                                          G4Molecule* molecule)
{
  MoleculeTable::iterator it = fMoleculeTable.find(name);

  G4int id = molecule->GetMoleculeID();
  MoleculeTablePerID::iterator it_perID = fMoleculeTablePerID.find(id);

  if (it == fMoleculeTable.end() && it_perID == fMoleculeTablePerID.end())
  {
    fMoleculeTable[name] = molecule;
    fMoleculeTablePerID[id] = molecule;
  }
  else if(it != fMoleculeTable.end())
  {
    // exception
    G4ExceptionDescription description;
    description << "The molecule model " << name
                << " was already recorded in the table" << G4endl;
    G4Exception("G4MoleculeTable::CreateMoleculeModel", "MODEL_ALREADY_CREATED",
                FatalException, description);
  }
  else if(it != fMoleculeTable.end())
  {
    // exception
    G4ExceptionDescription description;
    description << "The molecule model " << name
                << " was already recorded per ID (" << id<<") in the table"
                << G4endl;
    G4Exception("G4MoleculeTable::CreateMoleculeModel", "MODEL_ALREADY_CREATED",
                FatalException, description);
  }
}

G4Molecule* G4MoleculeTable::GetMoleculeModel(const G4String& name, bool mustExist)
{
  MoleculeTable::iterator it = fMoleculeTable.find(name);
  G4Molecule* molecule(0);
  if (it != fMoleculeTable.end())
  {
    molecule = it->second;
  }
  else if(mustExist)
  {
    // exception
    G4ExceptionDescription description;
    description << "The molecule model " << name
                << " was already recorded in the table" << G4endl;
    G4Exception("G4MoleculeTable::CreateMoleculeModel", "MODEL_ALREADY_CREATED",
                FatalException, description);
  }
  return molecule;
}

G4Molecule* G4MoleculeTable::GetMoleculeModel(G4int id)
{
  MoleculeTablePerID::iterator it = fMoleculeTablePerID.find(id);
  G4Molecule* molecule(0);
  if (it != fMoleculeTablePerID.end())
  {
    molecule = it->second;
  }
  else
  {
    // exception
    return 0;
  }
  return molecule;
}

void G4MoleculeTable::Insert(G4MoleculeDefinition* moleculeDefinition)
{

  const G4String& name = moleculeDefinition->GetName();
  MoleculeDefTable::iterator it = fMoleculeDefTable.find(name);
  if (it == fMoleculeDefTable.end())
  {
    fMoleculeDefTable[name] = moleculeDefinition;
  }
  else
  {
    // exception
    G4ExceptionDescription description;
    description << "The molecule definition " << name
                << " was already recorded in the table" << G4endl;
    G4Exception("G4MoleculeTable::CreateMoleculeDefinition",
                "DEFINITION_ALREADY_CREATED", FatalException, description);
  }
}
