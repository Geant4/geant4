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
// $Id: G4PolarizedPEEffectModel.hh 96114 2016-03-16 18:51:33Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4PolarizedPEEffectModel
//
// Author:        Andreas Schaelicke & Karim Laihem
//
// Creation date: 22.02.2007
//
// Modifications:
// 02.08.2007 : adapt to design change in version 4.9 (AS)
//
//
// Class Description:
//
// Implementation of polarization transfer in Photoelectric Effect
// 

// -------------------------------------------------------------------
//

#ifdef NOIONIZATIONAS
#define G4PolarizedPEEffectModel_h 1
#endif

#ifndef G4PolarizedPEEffectModel_h
#define G4PolarizedPEEffectModel_h 1

#include "G4PEEffectFluoModel.hh"
#include "G4StokesVector.hh"


class G4PolarizedPEEffectCrossSection;

class G4PolarizedPEEffectModel : public G4PEEffectFluoModel
{

public:

  explicit G4PolarizedPEEffectModel(const G4ParticleDefinition* p = nullptr,
			   const G4String& nam = "Polarized-PhotoElectric");

  void Initialise(const G4ParticleDefinition* pd, 
                  const G4DataVector& dv) override;
  virtual ~G4PolarizedPEEffectModel();

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;

  // polarized routines 
  /*
  inline void SetTargetPolarization(const G4ThreeVector & pTarget);
  inline void SetBeamPolarization(const G4ThreeVector & pBeam);
  inline const G4ThreeVector & GetTargetPolarization() const;
  inline const G4ThreeVector & GetBeamPolarization() const;
  inline const G4ThreeVector & GetFinalGammaPolarization() const;
  inline const G4ThreeVector & GetFinalElectronPolarization() const;
  */
private:

  // hide assignment operator
  G4PolarizedPEEffectModel & 
    operator=(const  G4PolarizedPEEffectModel &right) = delete;
  G4PolarizedPEEffectModel(const  G4PolarizedPEEffectModel&) = delete;

  G4PolarizedPEEffectCrossSection * crossSectionCalculator;
  // incomming
  G4StokesVector theBeamPolarization;    // photon
  // outgoing
  G4StokesVector finalElectronPolarization;

  G4int verboseLevel;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#endif
