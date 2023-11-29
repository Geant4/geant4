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

#include "G4MoleculeTableMessenger.hh"
#include "G4UIcmdWithAString.hh"
#include <G4SystemOfUnits.hh>
#include <G4UIcmdWithoutParameter.hh>
#include "G4MoleculeTable.hh"
#include "G4MoleculeDefinition.hh"
#include "G4MolecularConfiguration.hh"
#include "G4ParticleTable.hh"
//------------------------------------------------------------------------------

G4MoleculeTableMessenger::G4MoleculeTableMessenger()
  : G4UImessenger()
  , fpPrintTable(new G4UIcmdWithoutParameter("/chem/PrintSpeciesTable", this))
  , fpSpecies(new G4UIcmdWithAString("/chem/species", this))
{}

G4MoleculeTableMessenger::~G4MoleculeTableMessenger() = default;

//------------------------------------------------------------------------------
void G4MoleculeTableMessenger::SetNewValue(G4UIcommand* command,
                                           G4String newValue)
{
  if(command == fpPrintTable.get())
  {
    G4MolecularConfiguration::PrintAll();
  }
  if(command == fpSpecies.get())
  {
    std::istringstream iss(newValue);

    G4String speciesName;
    iss >> speciesName;

    G4String marker;
    iss >> marker;  // must be [

    if(marker != "[")
    {
      G4ExceptionDescription description;
      description << " marker : " << marker << G4endl;
      G4Exception("G4MoleculeTableMessenger::SetNewValue",
                  "FAIL_SPECIES_DEFINITION04", FatalException, description);
    }

    G4String speciesDefName;
    iss >> speciesDefName;

    iss >> marker;  // peut etre |

    G4int charge = 0;
    if(marker == "|")
    {
      iss >> charge;
    }
    iss >> marker;  // peut etre |
    G4double diffusion_coefficient = 0;
    if(marker == "|")
    {
      iss >> diffusion_coefficient;
    }
    iss >> marker;  // peut etre |
    G4double VanDerVaalsRadius = 0;
    if(marker == "|")
    {
      iss >> VanDerVaalsRadius;  // in nm
    }

    auto pConf =
      G4MolecularConfiguration::GetMolecularConfiguration(speciesName);
    if(pConf != nullptr)
    {
      pConf->UnFinalize();
      if(VanDerVaalsRadius != 0)
      {
        pConf->SetVanDerVaalsRadius(VanDerVaalsRadius * nm);
      }
      if(diffusion_coefficient != 0)
      {
        pConf->SetDiffusionCoefficient(diffusion_coefficient * (m2 / s));
      }
    }
    else
    {
      auto particleDef = dynamic_cast<G4MoleculeDefinition*>(
        G4ParticleTable::GetParticleTable()->FindParticle(speciesDefName));
      if(particleDef != nullptr)
      {
        auto molConf =
          G4MolecularConfiguration::GetOrCreateMolecularConfiguration(
            particleDef, charge);
        if(molConf == nullptr)
        {
          G4ExceptionDescription description;
          description << "This molecule has not been defined" << G4endl;
          G4Exception("G4MoleculeTableMessenger::SetNewValue",
                      "FAIL_SPECIES_DEFINITION02", FatalException, description);
        }
        molConf->UnFinalize();
        if(VanDerVaalsRadius != 0)
        {
          molConf->SetVanDerVaalsRadius(VanDerVaalsRadius * nm);
        }
        if(diffusion_coefficient != 0)
        {
          molConf->SetDiffusionCoefficient(diffusion_coefficient * (m2 / s));
        }

        auto usedName = molConf->GetUserID();
        if(usedName != "")
        {
          molConf->PrintState();
          G4ExceptionDescription description;
          description << "This molecule has been defined by the name : "
                      << usedName << " . Please, use this name." << G4endl;
          G4Exception("G4MoleculeTableMessenger::SetNewValue",
                      "FAIL_SPECIES_DEFINITION", FatalException, description);
        }
        else
        {
          molConf->SetUserID(speciesName);
        }
      }
      else
      {
        auto speciesDef =
          new G4MoleculeDefinition(speciesDefName, /*mass no needed*/ 0,
                                   /*D*/ diffusion_coefficient * (m * m / s),
                                   /*charge*/ charge, 1,
                                   /*electronL*/ 0,
                                   /*radius*/ VanDerVaalsRadius * nm);
        G4bool alreadyCreated(false);
        G4MolecularConfiguration::CreateMolecularConfiguration(
          speciesName, speciesDef, alreadyCreated);
      }
    }
  }
}
