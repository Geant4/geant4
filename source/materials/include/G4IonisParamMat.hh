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
// $Id: G4IonisParamMat.hh 106243 2017-09-26 01:56:43Z gcosmo $
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
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Threading.hh"
#include "G4DensityEffectCalc.hh"

class G4Material;                        // forward declaration
class G4DensityEffectData;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

class G4IonisParamMat  // with description
{
public:

  G4IonisParamMat(const G4Material*); 
  ~G4IonisParamMat();

  // parameters for mean energy loss calculation:
  inline
  G4double  GetMeanExcitationEnergy()   const {return fMeanExcitationEnergy;};

  void      SetMeanExcitationEnergy(G4double value);
  G4double  FindMeanExcitationEnergy(const G4Material*) const;

  inline
  G4double  GetLogMeanExcEnergy()       const {return fLogMeanExcEnergy;};
  inline
  G4double* GetShellCorrectionVector()  const {return fShellCorrectionVector;};
  inline
  G4double  GetTaul()                   const {return fTaul;};
    
  // parameters of the Sternheimer 3-part approximate density correction:
  inline
  G4double  GetPlasmaEnergy()           const {return fPlasmaEnergy;};
  inline
  G4double  GetAdjustmentFactor()       const {return fAdjustmentFactor;};
  inline
  G4double  GetCdensity()               const {return fCdensity;};
  inline
  G4double  GetMdensity()               const {return fMdensity;};
  inline
  G4double  GetAdensity()               const {return fAdensity;};
  inline
  G4double  GetX0density()              const {return fX0density;};
  inline
  G4double  GetX1density()              const {return fX1density;};
  inline
  G4double  GetD0density()              const {return fD0density;};

  // user defined density correction parameterisation.  These values are
  // only used if the more precise Sternheimer exact form cannot be
  // calculated.
  void SetDensityEffectParameters(G4double cd, G4double md, G4double ad,
                                  G4double x0, G4double x1, G4double d0);

  // defined density correction parameterisation via base material.
  // These values are only used if the more precise Sternheimer exact form
  // cannot be calculated.
  void SetDensityEffectParameters(const G4Material* bmat);

  /// Use the exact form of the density effect intead of the parameterization
  ///
  /// The "Sternheimer exact" form is calculated as described in
  ///
  /// R.M. Sternheimer et al. "Density Effect For The Ionization Loss
  /// of Charged Particles in Various Substances"
  /// Atom. Data Nucl. Data Tabl. 30 (1984) 261-271.
  ///
  /// It is not really the exact form, which requires detailed optical
  /// data to compute, but it is a better approximation than the popular
  /// parameterization introduced in Sternheimer's 1952 paper and used
  /// in the 1984 paper (among others).  If there were a significant
  /// number of substances for which the truly exact form could be
  /// calculated, I'd add another function SetExactDensityEffect().
  /// Alas, the list is only perhaps aluminum, water, and silicon.
  ///
  /// When this function is called, the atomic data for the material is
  /// used to set up the calculation. This involves finding a root of an
  /// equation and can fail. If it fails, this function returns false,
  /// a warning message will be printed, and we will fall back to the
  /// parameterization as though this function had not been called.
  ///
  /// Important for users:
  /// The effect is different for insulators and conductors. Materials are
  /// assumed to be non-conductors by default. To make a conductor:
  /// GetMaterialPropertiesTable()->AddConstProperty("conductor", 1)
  void SetSternheimerExactDensityEffect();

  // compute density correction as a function of the kinematic variable
  // x = log10(beta*gamma)  
  G4double DensityCorrection(G4double x);

  static G4DensityEffectData* GetDensityEffectData();

  // parameters of the energy loss fluctuation model:
  inline
  G4double  GetF1fluct()                const {return fF1fluct;};
  inline
  G4double  GetF2fluct()                const {return fF2fluct;};
  inline
  G4double  GetEnergy1fluct()           const {return fEnergy1fluct;};
  inline
  G4double  GetLogEnergy1fluct()        const {return fLogEnergy1fluct;};
  inline
  G4double  GetEnergy2fluct()           const {return fEnergy2fluct;};
  inline
  G4double  GetLogEnergy2fluct()        const {return fLogEnergy2fluct;};
  inline
  G4double  GetEnergy0fluct()           const {return fEnergy0fluct;};
  inline
  G4double  GetRateionexcfluct()        const {return fRateionexcfluct;};

  // parameters for ion corrections computations
  inline
  G4double  GetZeffective()             const {return fZeff;};
  inline
  G4double  GetFermiEnergy()            const {return fFermiEnergy;};
  inline
  G4double  GetLFactor()                const {return fLfactor;};
  inline
  G4double  GetInvA23()                 const {return fInvA23;};
    
  // parameters for Birks attenuation:
  inline
  void      SetBirksConstant(G4double value) {fBirks = value;}; 
  inline
  G4double  GetBirksConstant()         const {return fBirks;};

  // parameters for average energy per ion 
  inline
  void      SetMeanEnergyPerIonPair(G4double value) {fMeanEnergyPerIon = value;}; 
  inline
  G4double  GetMeanEnergyPerIonPair()         const {return fMeanEnergyPerIon;};
      
  G4IonisParamMat(__void__&);
  // Fake default constructor for usage restricted to direct object
  // persistency for clients requiring preallocation of memory for
  // persistifiable objects.

private:
    
  // Compute mean parameters : ExcitationEnergy,Shell corretion vector ...
  void ComputeMeanParameters();

  // Compute Sternheimer 3-part approximation parameters for the density
  // effect.  These are only used if the more precise Sternheimer exact form
  // cannot be used.
  void ComputeDensityEffect();

  // Compute parameters for the energy fluctuation model
  void ComputeFluctModel();

  // Compute parameters for ion parameterizations
  void ComputeIonParameters();

  // operators
  G4IonisParamMat& operator=(const G4IonisParamMat&) = delete;
  G4int operator==(const G4IonisParamMat&) const = delete;
  G4int operator!=(const G4IonisParamMat&) const = delete;
  G4IonisParamMat(const G4IonisParamMat&) = delete;

  //
  // data members
  //
  const G4Material* fMaterial;       // this material

  // parameters for mean energy loss calculation
  G4double  fMeanExcitationEnergy;          // 
  G4double  fLogMeanExcEnergy;              // 
  G4double* fShellCorrectionVector;         // shell correction coefficients
  G4double  fTaul;                          // lower limit of Bethe-Bloch formula

  // parameters of the Sternheimer approximate 3-part density correction
  G4double fCdensity;                      // mat.constant
  G4double fMdensity;                      // exponent
  G4double fAdensity;                      //
  G4double fX0density;                     //
  G4double fX1density;                     //
  G4double fD0density;

  // Stores atomic data necessary for computing the Sternheimer exact
  // density effect, as well as intermediate stages of the computation.
  G4DensityEffectCalcData * fCalcDensity;

  // If the computation of the Sternheimer exact density effect fails,
  // this flag is set to prevent wasting time on repeated failures.
  G4bool fCalcDensityFailure = false;

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
  G4double twoln10;
#ifdef G4MULTITHREADED
  static G4Mutex ionisMutex;
#endif
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
#endif
