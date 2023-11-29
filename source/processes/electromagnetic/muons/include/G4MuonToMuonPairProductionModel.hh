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
// File name:     G4MuonToMuonPairProductionModel
//
// Author:        Siddharth Yajaman on the base of Vladimir Ivantchenko code
//
// Creation date: 12.07.2022
//
// Modifications:
//

//
// Class Description:
//
// Implementation of mu+mu- pair production by muons
// R. P. Kokoulin, S. R. Kelner, and A. A. Petrukhin, 
// in Proceedings of the International Cosmic Ray Conference, 
// Salt Lake City, 1999, edited by D. Kieda et al. 
// (American Institute of Physics, Melville, NY, 2000), Vol. 2, p. 20
//
// It inherits from G4MuPairProductionModel
//
// -------------------------------------------------------------------
//

#ifndef G4MuonToMuonPairProductionModel_h
#define G4MuonToMuonPairProductionModel_h 1

#include "G4MuPairProductionModel.hh"
#include <vector>

class G4MuonToMuonPairProductionModel : public G4MuPairProductionModel
{
public:

  explicit G4MuonToMuonPairProductionModel(const G4ParticleDefinition* p = nullptr,
                                           const G4String& nam = "muToMuonPairProd");

  ~G4MuonToMuonPairProductionModel() = default;

  // hide assignment operator and copy constructor
  G4MuonToMuonPairProductionModel & operator=
  (const G4MuonToMuonPairProductionModel &right) = delete;
  G4MuonToMuonPairProductionModel(const G4MuonToMuonPairProductionModel&) = delete;

  G4double 
  ComputeDMicroscopicCrossSection(G4double tkin, G4double Z,
				  G4double pairEnergy) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*, 
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double maxEnergy) override;

  G4double U_func(G4double Z, G4double rho2, G4double xi, 
          G4double Y, G4double pairEnergy, const G4double B=183);

protected:

  void DataCorrupted(G4int Z, G4double logTkin) const override;

  G4ParticleDefinition* theMuonMinus = nullptr;
  G4ParticleDefinition* theMuonPlus = nullptr;

  G4double factorForCross;
  G4double minPairEnergy;
  G4double muonMass;
  G4double mueRatio;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
