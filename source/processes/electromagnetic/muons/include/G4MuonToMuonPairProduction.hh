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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4MuonToMuonPairProduction
//
// Author:        Siddharth Yajaman
//
// Creation date: 12.07.2022
//
// Modifications:

//
// Class Description:
//
// This class implements the muon pair production process for muons.
// it inherites from G4MuPairProduction.
//

// -------------------------------------------------------------------
//

#ifndef G4MuonToMuonPairProduction_h
#define G4MuonToMuonPairProduction_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MuPairProduction.hh"

class G4MuonToMuonPairProduction : public G4MuPairProduction
{
public:

  explicit G4MuonToMuonPairProduction(const G4String& processName = "muToMuonPairProd");

  ~G4MuonToMuonPairProduction() = default;

  G4MuonToMuonPairProduction & operator=
  (const G4MuonToMuonPairProduction &right) = delete;
  G4MuonToMuonPairProduction(const G4MuonToMuonPairProduction&) = delete;

  void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
			           const G4ParticleDefinition*) override;

  // print description in html
  void ProcessDescription(std::ostream&) const override;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
