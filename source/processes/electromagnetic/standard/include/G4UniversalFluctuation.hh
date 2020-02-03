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
// File name:     G4UniversalFluctuation
//
// Author:        V.Ivanchenko make a class with the Laszlo Urban model
//
// Creation date: 03.01.2002
//
// Modifications:
//
//
// Class Description:
//
// Implementation of energy loss fluctuations

// -------------------------------------------------------------------
//

#ifndef G4UniversalFluctuation_h
#define G4UniversalFluctuation_h 1


#include "G4VEmFluctuationModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4Poisson.hh"
#include <CLHEP/Random/RandomEngine.h>

class G4UniversalFluctuation : public G4VEmFluctuationModel
{

public:

  explicit G4UniversalFluctuation(const G4String& nam = "UniFluc");

  virtual ~G4UniversalFluctuation();

  virtual G4double SampleFluctuations(const G4MaterialCutsCouple*,
                                      const G4DynamicParticle*,
                                      G4double,
                                      G4double,
                                      G4double) override;

  virtual G4double Dispersion(const G4Material*,
                              const G4DynamicParticle*,
                              G4double,
                              G4double) override;

  virtual void InitialiseMe(const G4ParticleDefinition*) final;

  // Initialisation prestep
  virtual void SetParticleAndCharge(const G4ParticleDefinition*, 
                                    G4double q2) final;

private:

  inline void AddExcitation(CLHEP::HepRandomEngine* rndm, 
                            G4double a, G4double e, G4double& eav, 
                            G4double& eloss, G4double& esig2); 

  inline void SampleGauss(CLHEP::HepRandomEngine* rndm, 
                          G4double eav, G4double esig2, 
                          G4double& eloss); 

  // hide assignment operator
  G4UniversalFluctuation & operator=(const  G4UniversalFluctuation &right) = delete;
  G4UniversalFluctuation(const  G4UniversalFluctuation&) = delete;

  const G4ParticleDefinition* particle;
  const G4Material* lastMaterial;

  G4double particleMass;

  // Derived quantities
  G4double m_Inv_particleMass;
  G4double m_massrate;
  G4double chargeSquare;

  // data members to speed up the fluctuation calculation
  G4double ipotFluct;
  G4double electronDensity;
  
  G4double f1Fluct;
  G4double f2Fluct;
  G4double e1Fluct;
  G4double e2Fluct;
  G4double e1LogFluct;
  G4double e2LogFluct;
  G4double ipotLogFluct;
  G4double e0;
  G4double esmall;

  G4double e1,e2;

  G4double minNumberInteractionsBohr;
  G4double minLoss;
  G4double nmaxCont;
  G4double rate,a0,fw; 

  G4int     sizearray;
  G4double* rndmarray;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


inline void 
G4UniversalFluctuation::AddExcitation(CLHEP::HepRandomEngine* rndm, 
                                      G4double ax, G4double ex, G4double& eav,
                                      G4double& eloss, G4double& esig2) 
{
  if(ax > nmaxCont) {
    eav  += ax*ex;
    esig2 += ax*ex*ex;
  } else {
    G4int p = G4Poisson(ax);
    if(p > 0) { eloss += ((p + 1) - 2.*rndm->flat())*ex; }
  }
}

inline void 
G4UniversalFluctuation::SampleGauss(CLHEP::HepRandomEngine* rndm, 
                                    G4double eav, G4double esig2, 
                                    G4double& eloss)
{
  G4double x = eav;
  G4double sig = std::sqrt(esig2);
  if(eav < 0.25*sig) {
    x += (2.*rndm->flat() - 1.)*eav;
  } else {
    do { 
      x = G4RandGauss::shoot(rndm, eav, sig);
    } while (x < 0.0 || x > 2*eav);
    // Loop checking, 23-Feb-2016, Vladimir Ivanchenko
  }
  eloss += x;
} 

#endif

