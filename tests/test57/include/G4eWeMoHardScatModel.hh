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
// $Id: G4eWeMoHardScatModel.hh,v 1.1 2009-07-31 15:31:32 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eWeMoHardScatModel
//
// Author:        V. Grichine based on G4eCoulombScatteringModel 
//
// Creation date: 31.07.2009
//
//				       
//
// Class Description:
//
// Implementation of eCoulombScattering of pointlike charge particle 
// on Atomic Nucleus for interval of scattering anles in Lab system 
// thetaMin - ThetaMax, nucleus recoil is neglected.
//   The model based on analysis of J.M.Fernandez-Varea et al. 
// NIM B73(1993)447 originated from G.Wentzel Z.Phys. 40(1927)590 with 
// screening parameter from H.A.Bethe Phys. Rev. 89 (1953) 1256.
// 

// -------------------------------------------------------------------
//

#ifndef G4eWeMoHardScatModel_h
#define G4eWeMoHardScatModel_h 1

#include "G4VEmModel.hh"
#include "G4PhysicsTable.hh"
#include "globals.hh"
#include "G4NistManager.hh"

class G4ParticleChangeForGamma;
class G4ParticleDefinition;

class G4eWeMoHardScatModel : public G4VEmModel
{

public:

  G4eWeMoHardScatModel(const G4String& nam = "eCoulombScattering");
 
  virtual ~G4eWeMoHardScatModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
				G4double kinEnergy, 
				G4double Z, 
				G4double A, 
				G4double cut,
				G4double emax);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

  inline void SetRecoilThreshold(G4double eth);

protected:

  G4double CrossSectionPerAtom();

  G4double SampleCosineTheta();

  inline void DefineMaterial(const G4MaterialCutsCouple*);

  inline void SetupParticle(const G4ParticleDefinition*);

  inline void SetupKinematic(G4double kinEnergy, G4double cut);
  
  inline void SetupTarget(G4double Z, G4double kinEnergy); 

private:

  void ComputeMaxElectronScattering(G4double cut);

  // hide assignment operator
  G4eWeMoHardScatModel & operator=(const G4eWeMoHardScatModel &right);
  G4eWeMoHardScatModel(const  G4eWeMoHardScatModel&);

protected:
 
  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theElectron;
  const G4ParticleDefinition* thePositron;

  G4ParticleTable*          theParticleTable; 
  G4ParticleChangeForGamma* fParticleChange;
  G4NistManager*            fNistManager;
  const G4DataVector*       fCurrentCuts;

  const G4MaterialCutsCouple* fCurrentCouple;
  const G4Material*           fCurrentMaterial;
  const G4Element*            fCurrentElement;
  G4int                       fCurrentMaterialIndex;

  G4double                  fCoeff;
  G4double                  fCosThetaMin;
  G4double                  fCosThetaMax;
  G4double                  fCosTetMinNuc;
  G4double                  fCosTetMaxNuc;
  G4double                  fCosTetMaxNuc2;
  G4double                  fCosTetMaxElec;
  G4double                  fCosTetMaxElec2;
  G4double                  fq2Limit;
  G4double                  fRecoilThreshold;
  G4double                  fElecXSection;
  G4double                  fNucXSection;
  G4double                  feCut;

  // projectile

  const G4ParticleDefinition* fParticle;

  G4double                  fChargeSquare;
  G4double                  fSpin;
  G4double                  fMass;
  G4double                  fTkin;
  G4double                  fMom2;
  G4double                  fInvBeta2;
  //G4double                  kinFactor;
  G4double                  feTag;
  G4double                  fLowEnergyLimit;

  // target

  G4double                  fTargetZ;
  G4double                  fTargetMass;
  G4double                  fScreenZ;
  G4double                  fFormFactA;
  G4int                     fidxelm;
  G4int                     fiz;

private:

  G4double                  fAlpha2;
  G4double                  fFacLim;

  static G4double fScreenRSquare[100];
  static G4double fFormFactor[100];

  G4bool                    fInitialised;             
};

//////////////////////////////////////////////////////////////////////////////

inline
void G4eWeMoHardScatModel::DefineMaterial(const G4MaterialCutsCouple* cup) 
{ 
  if(cup != fCurrentCouple) 
  {
    fCurrentCouple        = cup;
    fCurrentMaterial      = cup->GetMaterial();
    fCurrentMaterialIndex = fCurrentCouple->GetIndex(); 
  }
}

//////////////////////////////////////////////////////////////////////////

inline 
void G4eWeMoHardScatModel::SetupParticle(const G4ParticleDefinition* p)
{
  // Initialise mass and charge

  if(p != fParticle) 
  {
    fParticle = p;
    fMass = fParticle->GetPDGMass();
    fSpin = fParticle->GetPDGSpin();
    G4double q = fParticle->GetPDGCharge()/eplus;
    fChargeSquare = q*q;
    fTkin = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////

inline void G4eWeMoHardScatModel::SetupKinematic( G4double ekin, 
						      G4double cut)
{
  if(ekin != fTkin || feCut != cut) 
  {
    fTkin = ekin;
    fMom2 = fTkin*(fTkin + 2.0*fMass);
    fInvBeta2 = 1.0 +  fMass*fMass/fMom2;
    fCosTetMinNuc = fCosThetaMin;
    fCosTetMaxNuc = fCosThetaMax;

    if(fMass < MeV && fCosThetaMin < 1.0 && ekin <= 10.*cut) 
    {
      fCosTetMinNuc = ekin*(fCosThetaMin + 1.0)/(10.*cut) - 1.0;
    }
    ComputeMaxElectronScattering(cut);
  }
}

////////////////////////////////////////////////////////////////////////
  
inline void G4eWeMoHardScatModel::SetupTarget(G4double Z, G4double e)
{
  if(Z != fTargetZ || e != feTag) 
  {
    feTag    = e; 
    fTargetZ = Z;
    fiz      = G4int(Z);

    if(fiz > 99) fiz = 99;
    fTargetMass   = fNistManager->GetAtomicMassAmu(fiz)*amu_c2;
    //  G4double m12 = mass*mass;
    //  G4double x   = 1.0 + mass/targetMass;
    //  kinFactor    = (1.0 +  m12/mom2)/mom2;

    fScreenZ = fScreenRSquare[fiz]/fMom2;

    //  if(iz > 1) {
    //    screenZ *=(1.13 + 3.76*Z*Z*invbeta2*alpha2);
    //    kinFactor /= (x*x);
    // }

    if(fiz > 1) fScreenZ *= (1.13 + 3.76*Z*Z*fInvBeta2*fAlpha2);

      //    screenZ *=(1.13 + std::min(0.5,3.76*Z*Z*invbeta2*alpha2));

    fFormFactA = fFormFactor[fiz]*fMom2;
    fCosTetMaxNuc2 = fCosTetMaxNuc;

    if(1 == fiz && fParticle == theProton && fCosTetMaxNuc2 < 0.0) 
    {
      fCosTetMaxNuc2 = 0.0;
    }
    fCosTetMaxElec2 = fCosTetMaxElec;
  } 
} 

////////////////////////////////////////////////////////////////////////////

inline void G4eWeMoHardScatModel::SetRecoilThreshold(G4double eth)
{
  fRecoilThreshold = eth;
}

//
//
////////////////////////////////////////////////////////////////////////////

#endif
