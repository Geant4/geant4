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
/// \file scavenger/include/ParserChemReaction.hh
/// \brief Definition of the scavenger::ParserChemReaction class

#ifndef SCAVENGER_ParserChemReaction_h
#define SCAVENGER_ParserChemReaction_h 1

#include "globals.hh"
#include <vector>
#include <array>

/// Parser to read user files defining chemical reactions and
/// scavengers (reaction with background)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace scavenger
{

class ParserChemReaction {
public:
  ParserChemReaction() = default;

  ~ParserChemReaction() = default;

  void ReadReactionFile(const G4String &fileName);

  inline auto GetListReactant1() { return fListReactant1; }

  inline auto GetListReactant2() { return fListReactant2; }

  inline auto GetListProduct() { return fListProduct; }

  inline auto GetListRate() { return fListRate; }

private:
  void AddReaction(const G4String &reactionString,
                   const G4String &type);

  void ReadReservoir(const G4String &reservoirString);

  G4double GetScavengerConcentration(const G4String &name);

  static void ReadReaction(const G4String &reactionString,
                           std::vector<G4String> &reactant,
                           std::vector<G4String> &product,
                           G4double &reactionRate);

  void ImplementReaction(const G4String &reactant1,
                         const G4String &reactant2,
                         const std::vector<G4String> &product,
                         const G4double &reactionRate,
                         const G4String &type);

  static void ReplaceString(G4String &aString, const G4String &from,
                            const G4String &to);

  std::array<std::vector<G4String>, 5> fListReactant1;
  std::array<std::vector<G4String>, 5> fListReactant2;
  std::array<std::vector<std::vector<G4String>>, 5> fListProduct;
  std::array<std::vector<G4double>, 5> fListRate;
  std::map<G4String, G4double> fReservoirConcentrationMap;
};

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
