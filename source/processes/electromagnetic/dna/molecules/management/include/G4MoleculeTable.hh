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
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4MOLECULETABLE_HH_
#define G4MOLECULETABLE_HH_

#include "globals.hh"
#include "G4MoleculeIterator.hh"
#include "G4ElectronOccupancy.hh"

#include <memory>

class G4MoleculeDefinition;
class G4MolecularConfiguration;
class G4MoleculeTableMessenger;

typedef G4MoleculeIterator<G4MoleculeDefinition> G4MoleculeDefinitionIterator;
typedef G4MoleculeIterator<G4MolecularConfiguration> G4ConfigurationIterator;

class G4MoleculeTable
{
public:
  static G4MoleculeTable* Instance();
  static G4MoleculeTable* GetMoleculeTable();
  virtual ~G4MoleculeTable();

  //____________________________________________________________________________
  // The methods below enable to create G4MoleculeDefinition &
  // G4MolecularConfiguration with a user identifier so that they can be retrieved
  // from this molecule table
  //

  //____________________________________________________________________________

  G4MoleculeDefinition* CreateMoleculeDefinition(const G4String& userIdentifier,
                                                 double diffusion_coefficient);

  //____________________________________________________________________________

  G4MolecularConfiguration*
  CreateConfiguration(const G4String& userIdentifier,
                      const G4MoleculeDefinition* molDef,
                      const G4String& configurationLabel,
                      const G4ElectronOccupancy& eOcc);

  G4MolecularConfiguration*
  CreateConfiguration(const G4String& userIdentifier,
                      G4MoleculeDefinition*,
                      int charge,
                      double diffusion_coefficient = -1);

  G4MolecularConfiguration*
  CreateConfiguration(const G4String& userIdentifier,
                      G4MoleculeDefinition*);

  G4MolecularConfiguration*
  CreateConfiguration(const G4String& userIdentifier,
                      G4MoleculeDefinition*,
                      const G4String& configurationLabel,
                      int charge = 0);

  //____________________________________________________________________________

  G4MoleculeDefinition* GetMoleculeDefinition(const G4String&,
                                              bool mustExist = true);

  G4MolecularConfiguration* GetConfiguration(const G4String&,
                                             bool mustExist = true);
  G4MolecularConfiguration* GetConfiguration(G4int id);

  //____________________________________________________________________________

  void Insert(G4MoleculeDefinition*);
  void Finalize(G4MoleculeDefinition*){}
  void Finalize();
  //____________________________________________________________________________

  G4MoleculeDefinitionIterator GetDefintionIterator()
  {
    return G4MoleculeDefinitionIterator(this->fMoleculeDefTable);
  }

  G4ConfigurationIterator GetConfigurationIterator();

  void PrepareMolecularConfiguration();

  int GetNumberOfDefinedSpecies();

protected:
  G4MoleculeTable();

  static G4MoleculeTable* fpgMoleculeTable;
  typedef std::map<G4String, G4MoleculeDefinition*> MoleculeDefTable;

  MoleculeDefTable fMoleculeDefTable;
  std::unique_ptr<G4MoleculeTableMessenger> fMoleculeDefTableMessenger;

};

#endif /* G4MOLECULETABLE_HH_ */
