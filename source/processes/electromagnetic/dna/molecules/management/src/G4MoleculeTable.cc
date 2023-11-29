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
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeDefinition.hh"
#include "G4MoleculeTableMessenger.hh"

G4MoleculeTable* G4MoleculeTable::fpgMoleculeTable(0);

//------------------------------------------------------------------------------

G4MoleculeTable::G4MoleculeTable()
  : fMoleculeDefTableMessenger(new G4MoleculeTableMessenger())
{
}

//------------------------------------------------------------------------------

G4MoleculeTable::~G4MoleculeTable()
{
}

//------------------------------------------------------------------------------

G4MoleculeTable* G4MoleculeTable::Instance()
{
  if (!fpgMoleculeTable) fpgMoleculeTable = new G4MoleculeTable;
  return fpgMoleculeTable;
}

//------------------------------------------------------------------------------

G4MoleculeTable* G4MoleculeTable::GetMoleculeTable()
{
  return Instance();
}

//------------------------------------------------------------------------------

G4MoleculeDefinition*
G4MoleculeTable::CreateMoleculeDefinition(const G4String& name,
                                          double diffusion_coefficient)
{
  return new G4MoleculeDefinition(name, -1 /* mass*/,
                                  diffusion_coefficient);
}

//------------------------------------------------------------------------------

G4MoleculeDefinition*
G4MoleculeTable::GetMoleculeDefinition(const G4String& name,
                                       bool mustExist)
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
    description << "The molecule definition " << name
                << " was NOT recorded in the table" << G4endl;
    G4Exception("G4MoleculeTable::CreateMoleculeModel",
                "MOLECULE_DEFINITION_NOT_CREATED",
                FatalException,
                description);
  }
  return definition;
}

//------------------------------------------------------------------------------

G4MolecularConfiguration*
G4MoleculeTable::GetConfiguration(const G4String& name, bool mustExist)
{
  G4MolecularConfiguration* species =
      G4MolecularConfiguration::GetMolecularConfiguration(name);

  if(species == 0 && mustExist)
  {
    // exception
    G4ExceptionDescription description;
    description << "The configuration " << name
                << " was not recorded in the table" << G4endl;
    G4Exception("G4MoleculeTable::GetConfiguration",
                "CONF_NOT_CREATED",
                FatalException,
                description);
  }

  return species;
}

//------------------------------------------------------------------------------

G4MolecularConfiguration*
G4MoleculeTable::GetConfiguration(G4int id)
{
  G4MolecularConfiguration* species =
        G4MolecularConfiguration::GetMolecularConfiguration(id);

  return species;
}

//------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------

void G4MoleculeTable::PrepareMolecularConfiguration()
{
  MoleculeDefTable::iterator it = fMoleculeDefTable.begin();

  for(; it != fMoleculeDefTable.end() ; ++it)
  {
    G4MolecularConfiguration::GetOrCreateMolecularConfiguration(it->second);
  }
}

//------------------------------------------------------------------------------

G4MolecularConfiguration*
G4MoleculeTable::CreateConfiguration(const G4String& userIdentifier,
                                     G4MoleculeDefinition* molDef)
{
  bool alreadyCreated(false);

  G4MolecularConfiguration* molConf =
      G4MolecularConfiguration::CreateMolecularConfiguration(userIdentifier,
                                                             molDef,
                                                             alreadyCreated);

  return molConf;
}

//------------------------------------------------------------------------------

G4MolecularConfiguration*
G4MoleculeTable::CreateConfiguration(const G4String& userIdentifier,
                                     G4MoleculeDefinition* molDef,
                                     const G4String& configurationLabel,
                                     int charge)
{
  bool alreadyCreated(false);

  G4MolecularConfiguration* molConf =
      G4MolecularConfiguration::CreateMolecularConfiguration(userIdentifier,
                                                             molDef,
                                                             charge,
                                                             configurationLabel,
                                                             alreadyCreated);

  return molConf;
}

//------------------------------------------------------------------------------

G4MolecularConfiguration*
G4MoleculeTable::CreateConfiguration(const G4String& userIdentifier,
                                     G4MoleculeDefinition* molDef,
                                     int charge,
                                     double diffusion_coefficient)
{
  bool alreadyCreated(false);

  G4MolecularConfiguration* molConf =
        G4MolecularConfiguration::CreateMolecularConfiguration(userIdentifier,
                                                               molDef,
                                                               charge,
                                                               userIdentifier,
                                                               alreadyCreated);

  if(diffusion_coefficient!=-1) // TODO
  {
    molConf->SetDiffusionCoefficient(diffusion_coefficient);
  }
  return molConf;
}

//------------------------------------------------------------------------------

G4MolecularConfiguration*
G4MoleculeTable::CreateConfiguration(const G4String& userIdentifier,
                                     const G4MoleculeDefinition* molDef,
                                     const G4String& configurationLabel,
                                     const G4ElectronOccupancy& eOcc)
{
  bool alreadyCreated(false);

  G4MolecularConfiguration* molConf =
      G4MolecularConfiguration::CreateMolecularConfiguration(userIdentifier,
                                                             molDef,
                                                             configurationLabel,
                                                             eOcc,
                                                             alreadyCreated);

  return molConf;
}

//------------------------------------------------------------------------------

void G4MoleculeTable::Finalize()
{
  G4MolecularConfiguration::FinalizeAll();
}

//------------------------------------------------------------------------------

G4ConfigurationIterator G4MoleculeTable::GetConfigurationIterator()
{
  return G4ConfigurationIterator(G4MolecularConfiguration::GetUserIDTable());
}

//------------------------------------------------------------------------------

int G4MoleculeTable::GetNumberOfDefinedSpecies()
{
  return G4MolecularConfiguration::GetNumberOfSpecies();
}
