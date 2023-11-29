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
#ifndef MOLECULAR_CHROMOSOME_MAPPER_HH
#define MOLECULAR_CHROMOSOME_MAPPER_HH

#include "globals.hh"
#include "G4ThreeVector.hh"

#include "ChromosomeFactory.hh"

#include <map>
#include <vector>

class VirtualChromosome;

class ChromosomeMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ChromosomeMapper
{
 public:
  ChromosomeMapper();

  virtual ~ChromosomeMapper();

  [[maybe_unused]] VirtualChromosome* GetChromosome(
    const G4String&) const;

  [[maybe_unused]] VirtualChromosome* GetChromosome(
    const G4ThreeVector& p) const;

  [[maybe_unused]] G4String GetCurrentChromosomeKey(
    const G4ThreeVector&) const;

  [[maybe_unused]] std::vector<G4String> GetChromosomeKeys() const;

  inline G4int GetNumberOfChromosomes()
  {
    return (G4int) fChromosomes.size();
  };

  void AddChromosome(const G4String&, const std::vector<G4String>&);

  void Test();

  void SavePlotData(const G4String& filename);

 private:
  std::map<uint32_t, VirtualChromosome*> fChromosomes;
  ChromosomeMessenger* fpChromosomeMessenger;
  ChromosomeFactory* fpChromosomeFactory;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_CHROMOSOME_MAPPER_HH
