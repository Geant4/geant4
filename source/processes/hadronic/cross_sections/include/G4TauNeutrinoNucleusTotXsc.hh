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
// 14.08.17 V. Grichine
//
// 31.10.22 V. Grichine - extension to tau-neutrinos based on mu_nu XS energy scaled

#ifndef G4TauNeutrinoNucleusTotXsc_h
#define G4TauNeutrinoNucleusTotXsc_h


#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"

class G4ParticleDefinition;

class G4TauNeutrinoNucleusTotXsc : public G4VCrossSectionDataSet
{
public:
   
  G4TauNeutrinoNucleusTotXsc();
  ~G4TauNeutrinoNucleusTotXsc();

  virtual
  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A, const G4Element*, const G4Material*);
  
  virtual
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int , const G4Material*){ return true; };


  // virtual G4double GetElementCrossSection(const G4DynamicParticle*, G4int Z, const G4Material*);
  virtual G4double GetElementCrossSection(const G4DynamicParticle* dynPart,
					       G4int Z,
					  const G4Material* mat);
  virtual
  G4double GetIsoCrossSection(const G4DynamicParticle* aPart, G4int Z, G4int A,  
			      const G4Isotope*,
			      const G4Element*,
			      const G4Material*);

  G4int GetEnergyIndex(G4double energy);
  G4double GetNuMuTotCsXsc(G4int index, G4double energy, G4int Z, G4int A);
  G4double GetANuMuTotCsXsc(G4int index, G4double energy, G4int Z, G4int A);

  G4double GetNuMuTotCsArray(G4int index);
  G4double GetANuMuTotCsArray(G4int index);

  void SetCutEnergy(G4double ec){fCutEnergy=ec;};
  G4double GetCutEnergy(){return fCutEnergy;};

  void SetBiasingFactor(G4double bf){fBiasingFactor=bf;};
  G4double GetBiasingFactor(){return fBiasingFactor;};

  G4double GetTotXsc(){return fTotXsc;};
  G4double GetCcTotRatio(){return fCcTotRatio;};
  G4double GetQEratio(){return fQEratio;};

protected:

  G4double fCofXsc;    // 2*Gf*Gf*MeC2/pi
  G4double fSin2tW;    // sin^2theta_Weinberg
  G4double fCofS, fCofL;
  G4double fCutEnergy; // minimal recoil electron energy detected
  G4double fBiasingFactor; // biasing xsc up
  G4double fTotXsc, fCcTotRatio, fCcFactor, fNcFactor, fQEratio;
  G4double fEmc, fEtc, fDtc;
  G4int fIndex;

  static const G4double fNuMuEnergy[50];
  static const G4double fNuMuInXsc[50];
  static const G4double fNuMuQeXsc[50];
  static const G4double fANuMuInXsc[50];
  static const G4double fANuMuQeXsc[50];

  // G4ParticleDefinition* theMuonMinus;
  // G4ParticleDefinition* theMuonPlus;
};

#endif
