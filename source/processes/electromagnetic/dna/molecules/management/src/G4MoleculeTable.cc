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
  MoleculeTable::iterator it = fMoleculeTable.find(name);
  G4Molecule* molecule(0);
  if (it == fMoleculeTable.end())
  {
    molecule = new G4Molecule(molDef, charge);
    if (diffusion_coefficient != -1) molecule->SetDiffusionCoefficient(
        diffusion_coefficient);
    fMoleculeTable[name] = molecule;
  }
  else
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

G4Molecule* G4MoleculeTable::CreateMoleculeModel(const G4String& name,
                                                 G4MoleculeDefinition* molDef)
{
  MoleculeTable::iterator it = fMoleculeTable.find(name);
  G4Molecule* molecule(0);
  if (it == fMoleculeTable.end())
  {
    molecule = new G4Molecule(molDef);
    fMoleculeTable[name] = molecule;
  }
  else
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
  if (it == fMoleculeTable.end())
  {
    fMoleculeTable[name] = molecule;
  }
  else
  {
    // exception
    G4ExceptionDescription description;
    description << "The molecule model " << name
                << " was already recorded in the table" << G4endl;
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
