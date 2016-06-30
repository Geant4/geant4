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
// $Id: G4PolarizedAnnihilationModel.hh 96114 2016-03-16 18:51:33Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4PolarizedAnnihilationModel
//
// Author:        Andreas Schaelicke and Pavel Starovoitov
//
// Creation date: 01.05.2005
//
// Modifications:
// 18-07-06 use newly calculated cross sections (P. Starovoitov)
// 21-08-06 update interface to geant4.8.1 (A. Schaelicke)
// 10-07-07 copied Initialise() method from G4eeToTwoGammaModel to provide a  
//          ParticleChangeForGamma object (A. Schaelicke)
//
//
// Class Description:
//
// Implementation of polarized gamma Annihilation scattering on free electron
// 

// -------------------------------------------------------------------
//

#ifndef G4PolarizedAnnihilationModel_h
#define G4PolarizedAnnihilationModel_h 1

#include "G4eeToTwoGammaModel.hh"
#include "G4StokesVector.hh"

class G4ParticleChangeForGamma;
class G4PolarizedAnnihilationCrossSection;

class G4PolarizedAnnihilationModel : public G4eeToTwoGammaModel
{

public:

  explicit G4PolarizedAnnihilationModel(const G4ParticleDefinition* p = nullptr, 
			const G4String& nam = "Polarized-Annihilation");

  virtual ~G4PolarizedAnnihilationModel();

  virtual void Initialise(const G4ParticleDefinition*, 
			  const G4DataVector&) override;
  virtual G4double ComputeCrossSectionPerElectron(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double cut,
                                      G4double emax) override;
  void ComputeAsymmetriesPerElectron(G4double gammaEnergy,
				     G4double & valueX,
				     G4double & valueA,
				     G4double & valueT);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;

  // polarized routines 
  inline void SetTargetPolarization(const G4ThreeVector & pTarget);
  inline void SetBeamPolarization(const G4ThreeVector & pBeam);
  inline const G4ThreeVector & GetTargetPolarization() const;
  inline const G4ThreeVector & GetBeamPolarization() const;
  inline const G4ThreeVector & GetFinalGamma1Polarization() const;
  inline const G4ThreeVector & GetFinalGamma2Polarization() const;
private:

  // hide assignment operator
  G4PolarizedAnnihilationModel & 
    operator=(const  G4PolarizedAnnihilationModel &right) = delete;
  G4PolarizedAnnihilationModel(const  G4PolarizedAnnihilationModel&) = delete;

  G4PolarizedAnnihilationCrossSection * crossSectionCalculator;
  // incomming
  G4StokesVector theBeamPolarization;    // positron
  G4StokesVector theTargetPolarization;  // electron
  // outgoing
  G4StokesVector finalGamma1Polarization;
  G4StokesVector finalGamma2Polarization;

  G4int verboseLevel;

  G4ParticleChangeForGamma* gParticleChange;
  G4bool gIsInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4PolarizedAnnihilationModel::SetTargetPolarization(const G4ThreeVector & pTarget)
{
  theTargetPolarization = pTarget;
}
inline void G4PolarizedAnnihilationModel::SetBeamPolarization(const G4ThreeVector & pBeam)
{
  theBeamPolarization = pBeam;
}
inline const G4ThreeVector & G4PolarizedAnnihilationModel::GetTargetPolarization() const
{
  return theTargetPolarization;
}
inline const G4ThreeVector & G4PolarizedAnnihilationModel::GetBeamPolarization() const
{
  return theBeamPolarization;
}
inline const G4ThreeVector &  G4PolarizedAnnihilationModel::GetFinalGamma1Polarization() const
{
  return finalGamma1Polarization;
}
inline const G4ThreeVector & G4PolarizedAnnihilationModel::GetFinalGamma2Polarization() const
{
  return finalGamma2Polarization;
}

#endif
