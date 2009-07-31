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
// $Id: G4WeMoSoftMscModel.hh,v 1.3 2009-07-31 15:31:32 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4WeMoSoftMscModel
//
// Author:        V.Grichine 
//
// Creation date: 31.07.2009 from G4WentzelVIModel
//
// Modifications:
//
//
// Class Description:
//
// Implementation of the model of multiple scattering based on
// G.Wentzel, Z. Phys. 40 (1927) 590.
// H.W.Lewis, Phys Rev 78 (1950) 526.
// J.M. Fernandez-Varea et al., NIM B73 (1993) 447.
// L.Urban, CERN-OPEN-2006-077.
// 
////////////////////////////////////////////////////////////////////////////////

#ifndef G4WeMoSoftMscModel_h
#define G4WeMoSoftMscModel_h 1

//////////////////////////////////////////////////////////////////////////////////

#include "G4VMscModel.hh"
#include "G4PhysicsTable.hh"
#include "G4MscStepLimitType.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4NistManager.hh"

class G4LossTableManager;
class G4ParticleChangeForMSC;
class G4ParticleDefinition;

///////////////////////////////////////////////////////////////////////////////////

class G4WeMoSoftMscModel : public G4VMscModel
{

public:

  G4WeMoSoftMscModel(const G4String& nam = "WeMoSoftMsc");

  virtual ~G4WeMoSoftMscModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
					      G4double KineticEnergy,
					      G4double AtomicNumber,
					      G4double AtomicWeight=0., 
					      G4double cut = DBL_MAX,
					      G4double emax= DBL_MAX);

  virtual void SampleScattering(const G4DynamicParticle*, G4double safety);

  virtual G4double ComputeTruePathLengthLimit(const G4Track& track,
					      G4PhysicsTable* theLambdaTable,
					      G4double currentMinimalStep);

  virtual G4double ComputeGeomPathLength(G4double truePathLength);

  virtual G4double ComputeTrueStepLength(G4double geomStepLength);

private:

  G4double ComputeTransportXSectionPerAtom();

  G4double ComputeXSectionPerVolume();

  void ComputeMaxElectronScattering(G4double cut);

  inline G4double GetLambda(G4double kinEnergy);

  inline void SetupParticle(const G4ParticleDefinition*);

  inline void SetupKinematic(G4double kinEnergy, G4double cut);
  
  inline void SetupTarget(G4double Z, G4double kinEnergy);

  inline void DefineMaterial(const G4MaterialCutsCouple*);

  //  hide assignment operator

  G4WeMoSoftMscModel & operator=(const  G4WeMoSoftMscModel &right);
  G4WeMoSoftMscModel(const  G4WeMoSoftMscModel&);

  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theElectron;
  const G4ParticleDefinition* thePositron;

  G4ParticleChangeForMSC*   fParticleChange;

  G4PhysicsTable*           theLambdaTable;
  G4PhysicsTable*           theLambda2Table;
  G4LossTableManager*       theManager;
  const G4DataVector*       fCurrentCuts;

  G4NistManager*            fNistManager;

  G4double fNumLimit;
  G4double ftLimitMinFix;
  G4double fInvSqrt12;

  // cash

  G4double fPreKinEnergy;
  G4double feCut;
  G4double fLambda1;
  G4double fTruePathLength;
  G4double fGeomPathLength;
  G4double fLambdaEff;
  G4double fCurrentRange; 

  G4double fPar1;
  G4double fPar2;
  G4double fPar3;

  G4double fTransportXsc;
  std::vector<G4double> fXsc;
  std::vector<G4double> fProb;
  G4int    fnElements;

  G4int    fnBins;
  G4int    fnWarnings;
  G4int    fnWarnLimit;

  G4int    fCurrentMaterialIndex;

  const G4MaterialCutsCouple* fCurrentCouple;
  const G4Material* fCurrentMaterial;

  // single scattering parameters

  G4double fCoeff;
  G4double fCosThetaMin;
  G4double fCosThetaMax;
  G4double fCosTetMaxNuc;
  G4double fCosTetMaxNuc2;
  G4double fCosTetMaxElec;
  G4double fCosTetMaxElec2;
  G4double fq2Limit;
  G4double fAlpha2;

  // projectile

  const G4ParticleDefinition* fParticle;

  G4double fChargeSquare;
  G4double fSpin;
  G4double fMass;
  G4double fTkin;
  G4double fMom2;
  G4double fInvBeta2;
  G4double fKinFactor;
  G4double feTag;
  G4double fLowEnergyLimit;

  // target

  G4double fTargetZ;
  G4double fTargetMass;
  G4double fScreenZ;
  G4double fFormFactA;
  G4int    fiz;

  static G4double fScreenRSquare[100];
  static G4double fFormFactor[100];

  // flags

  G4bool   fInitialized;
  G4bool   fInside;
};

//////////////////////////////////////////////////////////////////////////////


inline
void G4WeMoSoftMscModel::DefineMaterial(const G4MaterialCutsCouple* cup) 
{ 
  if( cup != fCurrentCouple ) 
  {
    fCurrentCouple        = cup;
    fCurrentMaterial      = cup->GetMaterial();
    fCurrentMaterialIndex = fCurrentCouple->GetIndex(); 
  }
}

//////////////////////////////////////////////////////////////////////////////

inline
G4double G4WeMoSoftMscModel::GetLambda(G4double e)
{
  G4double x;

  if(theLambdaTable) 
  {
    G4bool b;
    x = ((*theLambdaTable)[fCurrentMaterialIndex])->GetValue(e, b);
  } 
  else 
  {
    x = CrossSection(fCurrentCouple, fParticle, e, 
		     (*fCurrentCuts)[fCurrentMaterialIndex]);
  }
  if( x > DBL_MIN ) x = 1./x;
  else              x = DBL_MAX;
  return x;
}

////////////////////////////////////////////////////////////////////////////////

inline void G4WeMoSoftMscModel::SetupParticle(const G4ParticleDefinition* p)
{
  // Initialise mass and charge

  if( p != fParticle ) 
  {
    fParticle     = p;
    fMass         = fParticle->GetPDGMass();
    fSpin         = fParticle->GetPDGSpin();
    G4double q   = fParticle->GetPDGCharge()/eplus;
    fChargeSquare = q*q;
    fTkin         = 0.0;
  }
}

////////////////////////////////////////////////////////////////////////////////

inline void G4WeMoSoftMscModel::SetupKinematic(G4double ekin, G4double cut)
{
  if( ekin != fTkin || feCut != cut) 
  {
    fTkin     = ekin;
    fMom2     = fTkin*(fTkin + 2.0*fMass);
    fInvBeta2 = 1.0 +  fMass*fMass/fMom2;

    fCosTetMaxNuc = fCosThetaMax;

    if(fMass < MeV && ekin <= 10.*cut) 
    {
      fCosTetMaxNuc = ekin*(fCosThetaMax + 1.0)/(10.*cut) - 1.0;
    }
    ComputeMaxElectronScattering(cut);
  } 
}

////////////////////////////////////////////////////////////////////////////////
  
inline void G4WeMoSoftMscModel::SetupTarget(G4double Z, G4double e)
{
  if( Z != fTargetZ || e != feTag ) 
  {
    feTag    = e; 
    fTargetZ = Z;
    fiz      = G4int(Z);

    if(fiz > 99) fiz = 99;

    fTargetMass      = fNistManager->GetAtomicMassAmu(fiz)*amu_c2;
    fScreenZ         = fScreenRSquare[fiz]/fMom2;
    G4double gamma2 = 1./(1.-1./fInvBeta2);
    G4double mu_c2  = (fMass*fTargetMass)/(fMass+fTargetMass);
    fKinFactor       = fCoeff*fTargetZ*fChargeSquare*fInvBeta2*fInvBeta2/(gamma2*mu_c2*mu_c2);
    /*
    G4double m12  = fMass*fMass;
    G4double x    = 1.0 + fMass/fTargetMass;
    fKinFactor  = fCoeff*fTargetZ*fChargeSquare*(1.0 +  m12/fMom2)/fMom2;
    fScreenZ = fScreenRSquare[iz]/fMom2;
    if(fiz > 1) 
    {
      fScreenZ *=(1.13 + 3.76*Z*Z*fAlpha2);
      fKinFactor /= (x*x);
    }
    */
    if(fiz > 1) fScreenZ *=(1.13 + 3.76*Z*Z*fAlpha2);

    //if(iz > 1) fScreenZ *=(1.13 + std::min(0.5,3.76*Z*Z*fInvBeta2*fAlpha2));

    fFormFactA      = fFormFactor[fiz]*fMom2;
    fCosTetMaxNuc2  = fCosTetMaxNuc;
    fCosTetMaxElec2 = fCosTetMaxElec;
  } 
} 

//
//
///////////////////////////////////////////////////////////////////////////////////



#endif

