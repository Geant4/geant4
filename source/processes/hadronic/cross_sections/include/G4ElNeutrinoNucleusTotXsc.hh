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
// The total muon (anti)neutrino-nucleus cross sections in the
// simplified form as A-multiplied nu_mu-nucleon cross-sections
//
// 24.04.20 V. Grichine
//
// (nu_e,anti_nu_e)-nucleus xsc

#ifndef G4ElNeutrinoNucleusTotXsc_h
#define G4ElNeutrinoNucleusTotXsc_h


#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"

class G4ParticleDefinition;

class G4ElNeutrinoNucleusTotXsc : public G4VCrossSectionDataSet
{
public:
   
  G4ElNeutrinoNucleusTotXsc();
  ~G4ElNeutrinoNucleusTotXsc();

  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int , G4int , const G4Element*, const G4Material*) override { return true; };
  
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int , const G4Material*) override { return false; };


  // virtual G4double GetElementCrossSection(const G4DynamicParticle*, G4int Z, const G4Material*);
  // G4double GetElementCrossSection(const G4DynamicParticle* dynPart,   G4int Z,  const G4Material* mat) override;

  G4double GetIsoCrossSection(const G4DynamicParticle* aPart, G4int Z, G4int A,  
			      const G4Isotope*,
			      const G4Element*,
			      const G4Material*) override;

  G4int GetEnergyIndex(G4double energy);
  G4double GetNuElTotCsXsc(G4int index, G4double energy);
  G4double GetANuElTotCsXsc(G4int index, G4double energy);

  G4double GetNuElTotCsArray(G4int index);
  G4double GetANuElTotCsArray(G4int index);

  void SetCutEnergy(G4double ec){fCutEnergy = ec;};
  G4double GetCutEnergy(){return fCutEnergy;};

  void SetBiasingFactor(G4double bf){fBiasingFactor=bf;};
  G4double GetBiasingFactor(){return fBiasingFactor;};

  G4double GetTotXsc(){return fTotXsc;};
  G4double GetCcTotRatio(){return fCcTotRatio;};

protected:

  G4double fCofXsc;    // 2*Gf*Gf*MeC2/pi
  G4double fSin2tW;    // sin^2theta_Weinberg
  G4double fCofS, fCofL;
  G4double fCutEnergy; // minimal recoil electron energy detected
  G4double fBiasingFactor; // biasing xsc up
  G4double fTotXsc, fCcTotRatio, fCcFactor, fNcFactor;

  G4int fIndex;

  static const G4double fNuElEnergy[50];
  static const G4double fNuElTotXsc[50];
  static const G4double fANuElTotXsc[50];

  const G4ParticleDefinition* theElectron;
  const G4ParticleDefinition* thePositron;
};

#endif
