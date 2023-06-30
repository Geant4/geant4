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

#ifndef CHEMISTRYWORLD_HH
#define CHEMISTRYWORLD_HH

#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UImessenger.hh"
#include "G4VChemistryWorld.hh"
#include "G4DNABoundingBox.hh"
#include "G4SystemOfUnits.hh"

class ChemistryWorld : public G4VChemistryWorld, public G4UImessenger {
public:
  ChemistryWorld();

  ~ChemistryWorld() override = default;

  void ConstructChemistryBoundary() override;

  void ConstructChemistryComponents() override;

  void SetNewValue(G4UIcommand *, G4String) override;

private:
  std::unique_ptr<G4UIdirectory> fpChemWoldDir;
  std::unique_ptr<G4UIcmdWithADouble> fpAddpH;
  std::unique_ptr<G4UIcmdWithAString> fpAddScavengerName;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit>  fpTargetVolume;
  G4double fpH = 7;
  G4double fHalfBox = 1.6 * um;
};

#endif
