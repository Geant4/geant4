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
// $Id: G4IonisParamMat.hh,v 1.18 2010-05-10 10:44:39 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// class description
//
// The class contains few (physical) quantities related to the Ionisation
// process, for a material defined by its pointer G4Material*
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// 09-07-98: data moved from G4Material (mma)
// 09-03-01: copy constructor and assignement operator in public (mma)
// 28-10-02: add setMeanExcitationEnergy (V.Ivanchenko)
// 27-09-07: add computation of parameters for ions (V.Ivanchenko)
// 04-03-08: add fBirks constant (mma)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#ifndef G4IonisParamMat_HH
#define G4IonisParamMat_HH

#include "G4ios.hh"
#include "globals.hh"

class G4Material;                        // forward declaration
class G4DensityEffectData;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

class G4IonisParamMat  // with description
{
public:

  G4IonisParamMat(G4Material*); 
  virtual ~G4IonisParamMat();

  //
  // retrieval methods
  //
     
  // parameters for mean energy loss calculation:
  G4double  GetMeanExcitationEnergy()   const {return fMeanExcitationEnergy;};
  void      SetMeanExcitationEnergy(G4double value);
  G4double  FindMeanExcitationEnergy(const G4String& chFormula);
  G4double  GetLogMeanExcEnergy()       const {return fLogMeanExcEnergy;};
  G4double* GetShellCorrectionVector()  const {return fShellCorrectionVector;};
  G4double  GetTaul()                   const {return fTaul;};
    
  // parameters of the density correction:
  G4double  GetPlasmaEnergy()           const {return fPlasmaEnergy;};
  G4double  GetAdjustmentFactor()       const {return fAdjustmentFactor;};
  G4double  GetCdensity()               const {return fCdensity;};
  G4double  GetMdensity()               const {return fMdensity;};
  G4double  GetAdensity()               const {return fAdensity;};
  G4double  GetX0density()              const {return fX0density;};
  G4double  GetX1density()              const {return fX1density;};
  G4double  GetD0density()              const {return fD0density;};
    
  // compute density correction as a function of the kinematic variable
  // x = log10(beta*gamma)  
  inline G4double DensityCorrection(G4double x);
  static G4DensityEffectData* GetDensityEffectData();

  // parameters of the energy loss fluctuation model:
  G4double  GetF1fluct()                const {return fF1fluct;};
  G4double  GetF2fluct()                const {return fF2fluct;};
  G4double  GetEnergy1fluct()           const {return fEnergy1fluct;};
  G4double  GetLogEnergy1fluct()        const {return fLogEnergy1fluct;};
  G4double  GetEnergy2fluct()           const {return fEnergy2fluct;};
  G4double  GetLogEnergy2fluct()        const {return fLogEnergy2fluct;};
  G4double  GetEnergy0fluct()           const {return fEnergy0fluct;};
  G4double  GetRateionexcfluct()        const {return fRateionexcfluct;};

  // parameters for ion corrections computations
  G4double  GetZeffective()             const {return fZeff;};
  G4double  GetFermiEnergy()            const {return fFermiEnergy;};
  G4double  GetLFactor()                const {return fLfactor;};
  G4double  GetInvA23()                 const {return fInvA23;};
    
  // parameters for Birks attenuation:
  void      SetBirksConstant(G4double value) {fBirks = value;}; 
  G4double  GetBirksConstant()         const {return fBirks;};

  // parameters for average energy per ion 
  void      SetMeanEnergyPerIonPair(G4double value) {fMeanEnergyPerIon = value;}; 
  G4double  GetMeanEnergyPerIonPair()         const {return fMeanEnergyPerIon;};
      
public:  // without description

  G4IonisParamMat(const G4IonisParamMat&);
  const G4IonisParamMat& operator=(const G4IonisParamMat&);          
  G4int operator==(const G4IonisParamMat&) const;
  G4int operator!=(const G4IonisParamMat&) const;

  G4IonisParamMat(__void__&);
  // Fake default constructor for usage restricted to direct object
  // persistency for clients requiring preallocation of memory for
  // persistifiable objects.

private:
    
  // Compute mean parameters : ExcitationEnergy,Shell corretion vector ...
  void ComputeMeanParameters();

  // Compute parameters for the density effect
  void ComputeDensityEffect();

  // Compute parameters for the energy fluctuation model
  void ComputeFluctModel();

  // Compute parameters for ion parameterizations
  void ComputeIonParameters();

private:

//
// data members
//
  G4Material* fMaterial;                    // this material

  // parameters for mean energy loss calculation
  G4double  fMeanExcitationEnergy;          // 
  G4double  fLogMeanExcEnergy;              // 
  G4double* fShellCorrectionVector;         // shell correction coefficients
  G4double  fTaul;                          // lower limit of Bethe-Bloch formula

  // parameters of the density correction
  G4double fCdensity;                      // mat.constant
  G4double fMdensity;                      // exponent
  G4double fAdensity;                      //
  G4double fX0density;                     //
  G4double fX1density;                     //
  G4double fD0density;

  G4double fPlasmaEnergy;
  G4double fAdjustmentFactor;

  // parameters of the energy loss fluctuation model
  G4double fF1fluct;                       
  G4double fF2fluct;                       
  G4double fEnergy1fluct;
  G4double fLogEnergy1fluct;
  G4double fEnergy2fluct;
  G4double fLogEnergy2fluct;
  G4double fEnergy0fluct;
  G4double fRateionexcfluct;

  // parameters for ion corrections computations
  G4double fZeff;
  G4double fFermiEnergy;
  G4double fLfactor;
  G4double fInvA23;
    
  // parameter for Birks attenuation
  G4double fBirks;
  // average energy per ion pair
  G4double fMeanEnergyPerIon;

  // static data created only once
  static G4DensityEffectData* fDensityData;
};

  // x = log10(beta*gamma)  
inline G4double G4IonisParamMat::DensityCorrection(G4double x)
{
  static const G4double twoln10 = 2.*std::log(10.);
  G4double y = 0.0;
  if(x < fX0density) {
    if(fD0density > 0.0) { y = fD0density*std::pow(10.,2*(x - fX0density)); }
  } else if(x >= fX1density) { y = twoln10*x - fCdensity; }
  else {y = twoln10*x - fCdensity + fAdensity*std::pow(fX1density - x, fMdensity);}
  return y;
}

#endif
