//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4eCoulombScatteringModel.hh,v 1.1 2006-03-21 15:40:54 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eCoulombScatteringModel
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 19.02.2006
//
// Modifications:
//
// Class Description:
//
// Implementation of eCoulombScattering of pointlike charge particle 
// on Atomic Nucleus for interval of scattering anles in Lab system 
// thetaMin - ThetaMax, nucleus recoil is neglected
// 

// -------------------------------------------------------------------
//

#ifndef G4eCoulombScatteringModel_h
#define G4eCoulombScatteringModel_h 1

#include "G4VEmModel.hh"
#include "G4PhysicsTable.hh"

class G4ParticleChangeForGamma;

class G4eCoulombScatteringModel : public G4VEmModel
{

public:

  G4eCoulombScatteringModel(G4double thetaMin = 0.0, G4double thetaMax = pi,
			   G4bool build = true, G4double tlim = TeV*TeV,
			   const G4String& nam = "eCoulombScattering");
 
  virtual ~G4eCoulombScatteringModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A, 
                                      G4double cut,
                                      G4double emax);

  virtual std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy);

private:

  G4double CalculateCrossSectionPerAtom(const G4ParticleDefinition*, 
					G4double kinEnergy, 
					G4double Z);

  // hide assignment operator
  G4eCoulombScatteringModel & operator=(const G4eCoulombScatteringModel &right);
  G4eCoulombScatteringModel(const  G4eCoulombScatteringModel&);

  G4ParticleChangeForGamma* fParticleChange;
  G4PhysicsTable*           theCrossSectionTable; 

  G4double                  coeff;
  G4double                  a0;
  G4double                  cosThetaMin;
  G4double                  cosThetaMax;
  G4double                  lowMomentum;
  G4double                  highMomentum;
  G4double                  q2Limit;

  G4int                     nbins;
  G4int                     nmax;

  G4bool                    buildTable;             
  G4bool                    isInitialised;             
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4eCoulombScatteringModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition* p,
                                             G4double kinEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  G4double x;
  if(theCrossSectionTable) {
    G4bool b;
    G4double mass = p->GetPDGMass();
    G4double momentum2 = kinEnergy*(kinEnergy + 2.0*mass);
    G4double e     = kinEnergy + mass;
    G4double beta2 = momentum2/(e*e);
    x = (((*theCrossSectionTable)[G4int(Z)]))->GetValue(momentum2, b)
      / (momentum2*beta2);
  } else x = CalculateCrossSectionPerAtom(p, kinEnergy, Z);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
