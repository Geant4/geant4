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
// $Id: G4PolarizedMollerBhabhaModel.hh 96114 2016-03-16 18:51:33Z gcosmo $
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4PolarizedMollerBhabhaModel
//
// Author:        A.Schaelicke on base of Vladimir Ivanchenko code
//
// Creation date: 10.11.2005
//
// Modifications:
//
// 20-08-05, modified interface (A.Schaelicke)
//
//
// Class Description:
//
// Physics implementation of polarized Bhabha/Moller scattering
//
// -------------------------------------------------------------------
//

#ifndef G4PolarizedMollerBhabhaModel_h
#define G4PolarizedMollerBhabhaModel_h 1

#include "G4VEmModel.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4StokesVector.hh"

class G4VPolarizedCrossSection;

class G4PolarizedMollerBhabhaModel : public G4MollerBhabhaModel
{

public:

  explicit G4PolarizedMollerBhabhaModel(const G4ParticleDefinition* p = nullptr,
			     const G4String& nam = "PolarizedMollerBhabha");

  virtual ~G4PolarizedMollerBhabhaModel();

  virtual G4double ComputeCrossSectionPerElectron(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double cut,
                                      G4double emax) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;


  // polarization access routines (may go if using a modified ParticleChange)

  void SetTargetPolarization(const G4ThreeVector & pTarget)
  {
    theTargetPolarization = pTarget;
  }
  void SetBeamPolarization(const G4ThreeVector & pBeam)
  {
    theBeamPolarization = pBeam;
  }
  const G4StokesVector & GetTargetPolarization()
  {
    return theTargetPolarization;
  }
  const G4StokesVector & GetBeamPolarization()
  {
    return theBeamPolarization;
  }
  const G4StokesVector & GetFinalElectronPolarization()
  {
    return fElectronPolarization;
  }
  const G4StokesVector &  GetFinalPositronPolarization()
  {
    return fPositronPolarization;
  }
private:

  // copy constructor and hide assignment operator
  G4PolarizedMollerBhabhaModel(G4PolarizedMollerBhabhaModel &) = delete;
  G4PolarizedMollerBhabhaModel & 
    operator=(const G4PolarizedMollerBhabhaModel &right) = delete;

  G4StokesVector theBeamPolarization;
  G4StokesVector theTargetPolarization;

  G4VPolarizedCrossSection * crossSectionCalculator;

  G4StokesVector fPositronPolarization;
  G4StokesVector fElectronPolarization;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
