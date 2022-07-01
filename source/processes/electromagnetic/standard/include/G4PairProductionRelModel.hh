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
// File name:     G4PairProductionRelModel
//
// Author:        Andreas Schaelicke
//
// Creation date: 02.04.2009
//
// Modifications:
// 28-05-18 New version with improved screening function approximation, improved
//          LPM function approximation, efficiency, documentation and cleanup. 
//          Corrected call to selecting target atom in the final state sampling. 
//          (M. Novak)
//
// Class Description:
//
// Implementation of gamma conversion to e+e- in the field of a nucleus 
// relativistic approximation
// 

// -------------------------------------------------------------------
//

#ifndef G4PairProductionRelModel_h
#define G4PairProductionRelModel_h 1

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4VEmModel.hh"
#include "G4Log.hh"
#include <vector>

class G4ParticleChangeForGamma;
class G4Pow;

class G4PairProductionRelModel : public G4VEmModel
{

public:

  explicit G4PairProductionRelModel(const G4ParticleDefinition* p = nullptr, 
                                    const G4String& nam = "BetheHeitlerLPM");
 
  ~G4PairProductionRelModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  void InitialiseLocal(const G4ParticleDefinition*, 
		       G4VEmModel* masterModel) override;

  G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
				      G4double kinEnergy, 
				      G4double Z, 
				      G4double A=0., 
				      G4double cut=0.,
				      G4double emax=DBL_MAX) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double maxEnergy) override;

  void SetupForMaterial(const G4ParticleDefinition*,
			const G4Material*,G4double) override;

  inline void   SetLPMflag(G4bool val) { fIsUseLPMCorrection = val;  }
  inline G4bool LPMflag() const        { return fIsUseLPMCorrection; }

  G4PairProductionRelModel & operator=
  (const G4PairProductionRelModel &right) = delete;
  G4PairProductionRelModel(const  G4PairProductionRelModel&) = delete;

protected:

  // for evaluating screening related functions
  inline void ComputePhi12(const G4double delta, 
			   G4double &phi1, G4double &phi2);
  inline G4double ScreenFunction1(const G4double delta);
  inline G4double ScreenFunction2(const G4double delta);
  inline void ScreenFunction12(const G4double delta, 
			       G4double &f1, G4double &f2);
  // helper methods for cross-section computation under different approximations
  G4double ComputeParametrizedXSectionPerAtom(G4double gammaEnergy, G4double Z);
  G4double ComputeXSectionPerAtom(G4double gammaEnergy, G4double Z);
  G4double ComputeDXSectionPerAtom(G4double eplusEnergy, G4double gammaEnergy, 
                                   G4double Z);
  G4double ComputeRelDXSectionPerAtom(G4double eplusEnergy, 
				      G4double gammaEnergy, G4double Z);

private:

  // for creating some data structure per Z 
  void InitialiseElementData();
  struct ElementData {
    G4double  fLogZ13;
    G4double  fCoulomb;
    G4double  fLradEl;
    G4double  fDeltaFactor;
    G4double  fDeltaMaxLow;
    G4double  fDeltaMaxHigh;
    G4double  fEtaValue;
    G4double  fLPMVarS1Cond;
    G4double  fLPMILVarS1Cond;
  };
  // for precomputing comp. intensive parts of LPM suppression functions and 
  // using them at run-time
  void InitLPMFunctions();
  void ComputeLPMGsPhis(G4double &funcGS, G4double &funcPhiS, 
                        const G4double varShat); 
  void GetLPMFunctions(G4double &lpmGs, G4double &lpmPhis, const G4double sval);
  void ComputeLPMfunctions(G4double &fXiS, G4double &fGS, G4double &fPhiS, 
                           const G4double eps, const G4double egamma, 
                           const G4int izet);
  struct LPMFuncs {
    LPMFuncs() : fIsInitialized(false), fISDelta(100.), fSLimit(2.) {}
    G4bool                 fIsInitialized;
    G4double               fISDelta;
    G4double               fSLimit;
    std::vector<G4double>  fLPMFuncG;
    std::vector<G4double>  fLPMFuncPhi;
  };

protected:
  static const G4int                gMaxZet;
  //
  static const G4double             gLPMconstant;
  //
  static const G4double             gXGL[8]; 
  static const G4double             gWGL[8];
  static const G4double             gFelLowZet[8];
  static const G4double             gFinelLowZet[8];
  //
  static const G4double             gXSecFactor;
  static const G4double             gEgLPMActivation;
  //  
  static std::vector<ElementData*>  gElementData;  
  static LPMFuncs                   gLPMFuncs;
  // 
  G4bool                            fIsUseLPMCorrection;
  G4bool                            fIsUseCompleteScreening;
  //
  G4double                          fLPMEnergy;
  //
  G4double                          fParametrizedXSectionThreshold;
  G4double                          fCoulombCorrectionThreshold;
  //
  G4Pow*                            fG4Calc;
  G4ParticleDefinition*             fTheGamma;
  G4ParticleDefinition*             fTheElectron;
  G4ParticleDefinition*             fThePositron;
  G4ParticleChangeForGamma*         fParticleChange;
};
//
// Bethe screening functions for the elastic (coherent) scattering:
// Bethe's phi1, phi2 coherent screening functions were computed numerically 
// by using (the universal) atomic form factors computed based on the Thomas-
// Fermi model of the atom (using numerical solution of the Thomas-Fermi 
// screening function instead of Moliere's analytical approximation). The 
// numerical results can be well approximated (better than Butcher & Messel 
// especially near the delta=1 limit) by:
// ## if delta <= 1.4 
//  phi1(delta) = 20.806 - delta*(3.190 - 0.5710*delta)   
//  phi2(delta) = 20.234 - delta*(2.126 - 0.0903*delta)
// ## if delta  > 1.4
//  phi1(delta) = phi2(delta) = 21.0190 - 4.145*ln(delta + 0.958)
// with delta = 136mc^2kZ^{-1/3}/[E(Eg-E)] = 136Z^{-1/3}eps0/[eps(1-eps)] where 
// Eg is the initial photon energy, E is the total energy transferred to one of 
// the e-/e+ pair, eps0 = mc^2/Eg and eps = E/Eg.

inline void G4PairProductionRelModel::ComputePhi12(const G4double delta,
						   G4double &phi1, 
						   G4double &phi2)
{
    if (delta > 1.4) {
      phi1 = 21.0190 - 4.145*G4Log(delta + 0.958);
      phi2 = phi1;
    } else {
      phi1 = 20.806 - delta*(3.190 - 0.5710*delta);
      phi2 = 20.234 - delta*(2.126 - 0.0903*delta);
    }
}

// Compute the value of the screening function 3*PHI1(delta) - PHI2(delta):
inline G4double G4PairProductionRelModel::ScreenFunction1(const G4double delta)
{
  return (delta > 1.4) ? 42.038 - 8.29*G4Log(delta + 0.958) 
                       : 42.184 - delta*(7.444 - 1.623*delta);
}

// Compute the value of the screening function 1.5*PHI1(delta) +0.5*PHI2(delta):
inline G4double G4PairProductionRelModel::ScreenFunction2(const G4double delta)
{
  return (delta > 1.4) ? 42.038 - 8.29*G4Log(delta + 0.958)
                       : 41.326 - delta*(5.848 - 0.902*delta);
}

// Same as ScreenFunction1 and ScreenFunction2 but computes them at once
inline void G4PairProductionRelModel::ScreenFunction12(const G4double delta, 
                                                     G4double &f1, G4double &f2)
{
  if (delta > 1.4) {
    f1 = 42.038 - 8.29*G4Log(delta + 0.958);
    f2 = f1;
  } else {
    f1 = 42.184 - delta*(7.444 - 1.623*delta);
    f2 = 41.326 - delta*(5.848 - 0.902*delta); 
  }
}

#endif
