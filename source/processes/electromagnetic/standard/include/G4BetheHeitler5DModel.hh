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
// File name:     G4BetheHeitler5DModel
//
// Authors:
// Igor Semeniouk and Denis Bernard,
// LLR, Ecole polytechnique & CNRS/IN2P3, 91128 Palaiseau, France
//
// Modifications:
// 27-10-17 New class (IgS)
// 19-01-18 version that calculates the pdf in the same way as in the fortran 
//          version (Denis Bernard)
// 04-06-18 Performance optimization of the final state sampling (M. Novak) 
//
// Class Description:
//
// Implementation of gamma convertion to e+e- in the field of a nucleus
//

// -------------------------------------------------------------------
//

#ifndef G4BetheHeitler5DModel_h
#define G4BetheHeitler5DModel_h 1

#include "G4BetheHeitlerModel.hh"

class G4IonTable;

class G4BetheHeitler5DModel : public G4BetheHeitlerModel
{

public:

  explicit G4BetheHeitler5DModel(const G4ParticleDefinition* p = nullptr,
                                 const G4String& nam = "BetheHeitler5D");

  virtual ~G4BetheHeitler5DModel();

  virtual void Initialise(const G4ParticleDefinition*,
			  const G4DataVector&) final;

  void SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                         const G4MaterialCutsCouple* couple,
                         const G4DynamicParticle* aDynamicGamma,
                         G4double, G4double) final;

  inline void SetVerbose(G4int val) { fVerbose = val; }

private:

  // hide assignment operator
  G4BetheHeitler5DModel& operator=(const G4BetheHeitler5DModel& right) = delete;
  G4BetheHeitler5DModel(const  G4BetheHeitler5DModel&) = delete;

  void BoostG4LorentzVector(const G4LorentzVector& p, const G4LorentzVector& q, 
                            G4LorentzVector& res) const;

  void BoostG4LorentzVector(const G4LorentzVector& p, const G4double qz,
                            const G4double qt, const G4double lffac,
                            const G4double imass, G4LorentzVector& res) const;

  G4double MaxDiffCrossSection(const G4double* par, G4double eZ,
                               G4double e, G4double loge) const;

  G4IonTable* theIonTable;

  G4int  fVerbose;
  G4int  fConversionType;
  G4bool iraw;
};
#endif
