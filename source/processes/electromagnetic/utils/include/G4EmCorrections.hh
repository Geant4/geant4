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
// $Id: G4EmCorrections.hh 103954 2017-05-04 11:29:22Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmCorrections
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 13.01.2005
//
// Modifications:
// 28.04.2006 General cleanup, add finite size corrections (V.Ivanchenko)
// 13.05.2006 Add corrections for ion stopping (V.Ivanhcenko)
// 20.05.2008 Removed Finite Size correction (V.Ivanchenko)
// 12.09.2008 Added inlined interfaces to effective charge (V.Ivanchenko)
// 19.04.2012 Fix reproducibility problem (A.Ribon)
//
// Class Description:
//
// This class provides calculation of EM corrections to ionisation
//

// -------------------------------------------------------------------
//

#ifndef G4EmCorrections_h
#define G4EmCorrections_h 1

#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4ionEffectiveCharge.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"

class G4VEmModel;
class G4PhysicsVector;
class G4IonTable;
class G4MaterialCutsCouple;
class G4LPhysicsFreeVector;
class G4Pow;

class G4EmCorrections
{

public:

  explicit G4EmCorrections(G4int verb);

  virtual ~G4EmCorrections();

  G4double HighOrderCorrections(const G4ParticleDefinition*,
                                const G4Material*,
                                G4double kineticEnergy,
                                G4double cutEnergy);

  G4double IonHighOrderCorrections(const G4ParticleDefinition*,
                                   const G4MaterialCutsCouple*,
                                   G4double kineticEnergy);

  G4double ComputeIonCorrections(const G4ParticleDefinition*,
                                 const G4Material*,
                                 G4double kineticEnergy);

  G4double IonBarkasCorrection(const G4ParticleDefinition*,
                               const G4Material*,
                               G4double kineticEnergy);

  G4double Bethe(const G4ParticleDefinition*,
                 const G4Material*,
                 G4double kineticEnergy);

  G4double SpinCorrection(const G4ParticleDefinition*,
                          const G4Material*,
                          G4double kineticEnergy);

  G4double KShellCorrection(const G4ParticleDefinition*,
                            const G4Material*,
                            G4double kineticEnergy);

  G4double LShellCorrection(const G4ParticleDefinition*,
                            const G4Material*,
                            G4double kineticEnergy);

  G4double ShellCorrection(const G4ParticleDefinition*,
                           const G4Material*,
                           G4double kineticEnergy);

  G4double ShellCorrectionSTD(const G4ParticleDefinition*,
                              const G4Material*,
                              G4double kineticEnergy);

  G4double DensityCorrection(const G4ParticleDefinition*,
                             const G4Material*,
                             G4double kineticEnergy);

  G4double BarkasCorrection(const G4ParticleDefinition*,
                            const G4Material*,
                            G4double kineticEnergy);

  G4double BlochCorrection(const G4ParticleDefinition*,
                           const G4Material*,
                           G4double kineticEnergy);

  G4double MottCorrection(const G4ParticleDefinition*,
                          const G4Material*,
                          G4double kineticEnergy);

  void AddStoppingData(G4int Z, G4int A, const G4String& materialName,
                       G4PhysicsVector* dVector);

  void InitialiseForNewRun();

  // effective charge correction using stopping power data
  G4double EffectiveChargeCorrection(const G4ParticleDefinition*,
                                     const G4Material*,
                                     G4double kineticEnergy);

  // effective charge of an ion
  inline G4double GetParticleCharge(const G4ParticleDefinition*,
                                    const G4Material*,
                                    G4double kineticEnergy);

  inline
  G4double EffectiveChargeSquareRatio(const G4ParticleDefinition*,
                                      const G4Material*,
                                      G4double kineticEnergy);

  // ionisation models for ions
  inline void SetIonisationModels(G4VEmModel* m1 = nullptr, 
                                  G4VEmModel* m2 = nullptr);

  inline G4int GetNumberOfStoppingVectors() const;

  inline void SetVerbose(G4int verb);

private:

  void Initialise();

  void BuildCorrectionVector();

  void SetupKinematics(const G4ParticleDefinition*,
                       const G4Material*,
                       G4double kineticEnergy);

  G4double KShell(G4double theta, G4double eta);

  G4double LShell(G4double theta, G4double eta);

  G4int Index(G4double x, const G4double* y, G4int n) const;

  G4double Value(G4double xv, G4double x1, G4double x2, 
                 G4double y1, G4double y2) const;

  G4double Value2(G4double xv, G4double yv, G4double x1, G4double x2,
                  G4double y1, G4double y2, G4double z11, G4double z21, 
		  G4double z12, G4double z22) const;

  // hide assignment operator
  G4EmCorrections & operator=(const G4EmCorrections &right) = delete;
  G4EmCorrections(const G4EmCorrections&) = delete;

  G4Pow* g4calc;

  static const G4double inveplus;
  static const G4double ZD[11];
  static const G4double UK[20];
  static const G4double VK[20];
  static G4double ZK[20];
  static const G4double Eta[29];
  static G4double CK[20][29];
  static G4double CL[26][28];
  static const G4double UL[26];
  static G4double VL[26];

  static G4LPhysicsFreeVector* BarkasCorr;
  static G4LPhysicsFreeVector* ThetaK;
  static G4LPhysicsFreeVector* ThetaL;

  G4double     alpha2;

  std::vector<const G4Material*> currmat;
  std::map< G4int, std::vector<G4double> > thcorr;
  size_t       ncouples;

  const G4ParticleDefinition* particle;
  const G4ParticleDefinition* curParticle;
  const G4Material*           material;
  const G4Material*           curMaterial;
  const G4ElementVector*      theElementVector;
  const G4double*             atomDensity;

  G4PhysicsVector*            curVector;

  G4IonTable*  ionTable;
  G4VEmModel*  ionLEModel;
  G4VEmModel*  ionHEModel;

  G4double  kinEnergy;
  G4double  mass;
  G4double  massFactor;
  G4double  eth;
  G4double  tau;
  G4double  gamma;
  G4double  bg2;
  G4double  beta2;
  G4double  beta;
  G4double  ba2;
  G4double  tmax;
  G4double  charge;
  G4double  q2;
  G4double  eCorrMin;
  G4double  eCorrMax;

  G4int     verbose;

  G4int     nK;
  G4int     nL;
  G4int     nEtaK;
  G4int     nEtaL;

  G4int     nbinCorr;
  G4int     numberOfElements;

  // Ion stopping data
  G4int                       nIons;
  G4int                       idx;
  G4int                       currentZ;
  std::vector<G4int>          Zion;
  std::vector<G4int>          Aion;
  std::vector<G4String>       materialName;

  std::vector<const G4ParticleDefinition*> ionList;

  std::vector<const G4Material*> materialList;
  std::vector<G4PhysicsVector*>  stopData;

  G4bool    isMaster;
  G4ionEffectiveCharge  effCharge;
};

inline G4int 
G4EmCorrections::Index(G4double x, const G4double* y, G4int n) const
{
  G4int iddd = n-1;
  // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  do {--iddd;} while (iddd>0 && x<y[iddd]);
  return iddd;
}

inline G4double G4EmCorrections::Value(G4double xv, G4double x1, G4double x2,
                                       G4double y1, G4double y2) const
{
  return y1 + (y2 - y1)*(xv - x1)/(x2 - x1);
}

inline G4double G4EmCorrections::Value2(G4double xv, G4double yv, 
                                        G4double x1, G4double x2,
                                        G4double y1, G4double y2,
                                        G4double z11, G4double z21, 
                                        G4double z12, G4double z22) const
{
  return (z11*(x2-xv)*(y2-yv) + z22*(xv-x1)*(yv-y1) +
          0.5*(z12*((x2-xv)*(yv-y1)+(xv-x1)*(y2-yv))+
               z21*((xv-x1)*(y2-yv)+(yv-y1)*(x2-xv))))
         / ((x2-x1)*(y2-y1));
}

inline void 
G4EmCorrections::SetIonisationModels(G4VEmModel* mod1, G4VEmModel* mod2)
{
  if(mod1) { ionLEModel = mod1; }
  if(mod2) { ionHEModel = mod2; }
}

inline G4int G4EmCorrections::GetNumberOfStoppingVectors() const
{
  return nIons;
}

inline G4double 
G4EmCorrections::GetParticleCharge(const G4ParticleDefinition* p,
                                   const G4Material* mat,
                                   G4double kineticEnergy)
{
  return effCharge.EffectiveCharge(p,mat,kineticEnergy);
}

inline G4double 
G4EmCorrections::EffectiveChargeSquareRatio(const G4ParticleDefinition* p,
                                            const G4Material* mat,
                                            G4double kineticEnergy)
{
  return effCharge.EffectiveChargeSquareRatio(p,mat,kineticEnergy);
}

inline void G4EmCorrections::SetupKinematics(const G4ParticleDefinition* p,
                                             const G4Material* mat,
                                             G4double kineticEnergy)
{
  if(kineticEnergy != kinEnergy || p != particle) {
    particle = p;
    kinEnergy = kineticEnergy;
    mass  = p->GetPDGMass();
    tau   = kineticEnergy / mass;
    gamma = 1.0 + tau;
    bg2   = tau * (tau+2.0);
    beta2 = bg2/(gamma*gamma);
    beta  = std::sqrt(beta2);
    ba2   = beta2/alpha2;
    G4double ratio = CLHEP::electron_mass_c2/mass;
    tmax  = 2.0*CLHEP::electron_mass_c2*bg2 
      /(1. + 2.0*gamma*ratio + ratio*ratio);
    charge  = p->GetPDGCharge()*inveplus;
    if(charge > 1.5) { charge = effCharge.EffectiveCharge(p,mat,kinEnergy); }
    q2 = charge*charge;
  }
  if(mat != material) {
    material = mat;
    theElementVector = material->GetElementVector();
    atomDensity  = material->GetAtomicNumDensityVector(); 
    numberOfElements = material->GetNumberOfElements();
  }
}

inline void G4EmCorrections::SetVerbose(G4int verb)
{
  verbose = verb;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
