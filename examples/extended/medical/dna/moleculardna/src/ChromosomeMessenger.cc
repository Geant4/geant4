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
#include "ChromosomeMessenger.hh"
#include "ChromosomeMapper.hh"
#include "UtilityFunctions.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChromosomeMessenger::ChromosomeMessenger(ChromosomeMapper* mapper)
  : G4UImessenger()
  , fpChromosomeMapper(mapper)
  , fpChromosomeDirectory(new G4UIdirectory("/chromosome/"))
  , fpAddChromosome(new G4UIcmdWithAString("/chromosome/add", this))
  , fpSavePlotData(new G4UIcmdWithAString("/chromosome/plotData", this))
{
  // Chromosomes
  fpChromosomeDirectory->SetGuidance("Commands for chromosome geometry.");
  fpAddChromosome->SetGuidance("Add a chromosomal region");
  fpAddChromosome->SetGuidance("format: shape name geometry-commands");
  fpAddChromosome->SetGuidance("shape: sphere || cyl");
  fpAddChromosome->SetGuidance("geometry-commands:");
  fpAddChromosome->SetGuidance("sphere: rad x y z unit [rx ry rz]");
  fpAddChromosome->SetGuidance("cyl: rad height x y z unit [rx ry rz]");
  fpAddChromosome->SetGuidance("Rotations are optional and in degrees");
  fpSavePlotData->SetGuidance("Save plot data to specified file");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChromosomeMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fpAddChromosome.get())
  {
    std::vector<G4String> values = utility::Split(newValue, ' ');
    G4String key                 = values[0];
    std::vector<G4String> commands;
    for(auto it = (values.begin() + 1); it != values.end(); it++)
    {
      commands.emplace_back(*it);
    }
    fpChromosomeMapper->AddChromosome(key, commands);
  }
  else if(command == fpSavePlotData.get())
  {
    fpChromosomeMapper->SavePlotData(newValue);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
