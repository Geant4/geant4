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
// $Id: G4eCoulombScatteringModel.hh,v 1.12 2007-08-14 17:10:33 vnivanch Exp $
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
// 01.08.06 V.Ivanchenko extend upper limit of table to TeV and review the
//          logic of building - only elements from G4ElementTable
// 08.08.06 V.Ivanchenko build internal table in ekin scale, introduce faclim
// 19.08.06 V.Ivanchenko add inline function ScreeningParameter and
//                       make some members protected
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

#ifndef G4eCoulombScatteringModel_h
#define G4eCoulombScatteringModel_h 1

#include "G4VEmModel.hh"
#include "G4PhysicsTable.hh"
#include "globals.hh"

class G4ParticleChangeForGamma;
class G4ParticleDefinition;

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

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

protected:

  virtual G4double CalculateCrossSectionPerAtom(
                                         const G4ParticleDefinition*, 
					 G4double kinEnergy, 
					 G4double Z, G4double A);

  inline void SetupKinematic(const G4ParticleDefinition*, G4double kinEnergy);
  
  inline void SetupTarget(G4double Z, G4double A, G4double kinEnergy); 

private:

  // hide assignment operator
  G4eCoulombScatteringModel & operator=(const G4eCoulombScatteringModel &right);
  G4eCoulombScatteringModel(const  G4eCoulombScatteringModel&);

protected:

  const G4ParticleDefinition* theProton;
  G4ParticleChangeForGamma* fParticleChange;

  G4double                  coeff;
  G4double                  constn;
  G4double                  cosThetaMin;
  G4double                  cosThetaMax;
  G4double                  cosTetMaxNuc;
  G4double                  q2Limit;

  // projectile
  const G4ParticleDefinition* particle;

  G4double                  chargeSquare;
  G4double                  mass;
  G4double                  tkin;
  G4double                  mom2;
  G4double                  invbeta2;

  // target
  G4double                  kineticEnergy;
  G4double                  targetZ;
  G4double                  targetA;
  G4double                  screenZ;
  G4double                  formfactA;

private:

  const G4ParticleDefinition* theElectron;
  const G4ParticleDefinition* thePositron;

  G4PhysicsTable*           theCrossSectionTable; 

  G4double                  a0;
  G4double                  lowKEnergy;
  G4double                  highKEnergy;
  G4double                  alpha2;
  G4double                  faclim;

  G4int                     nbins;
  G4int                     nmax;
  G4int                     index[100];

  G4bool                    buildTable;             
  G4bool                    isInitialised;             
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4eCoulombScatteringModel::ComputeCrossSectionPerAtom(
                const G4ParticleDefinition* p,
		G4double kinEnergy,
		G4double Z, G4double A,
		G4double, G4double)
{
  G4bool b;
  G4double x = 0.0;
  if(theCrossSectionTable) {
    x = std::exp((((*theCrossSectionTable)[index[G4int(Z)]]))
		 ->GetValue(kinEnergy, b));
  } else x = CalculateCrossSectionPerAtom(p, kinEnergy, Z, A);

  //  G4cout << "G4eCoulombScatteringModel:ComputeCSPerAtom e= " << kinEnergy 
  //	 << " Z= " << Z
  //       << "  CS= " << x << G4endl;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4eCoulombScatteringModel::SetupKinematic(
		const G4ParticleDefinition* p, 
		G4double ekin)
{
  if(ekin != tkin || p != particle) {
    particle = p;
    mass = particle->GetPDGMass();
    G4double q = particle->GetPDGCharge()/eplus;
    chargeSquare = q*q;
    tkin  = ekin;
    mom2  = tkin*(tkin + 2.0*mass);
    invbeta2 = 1.0 +  mass*mass/mom2;
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline void G4eCoulombScatteringModel::SetupTarget(G4double Z, G4double A, 
						   G4double e)
{
  if(Z != targetZ || A != targetA || e != kineticEnergy) {
    targetZ = Z;
    targetA = A;
    kineticEnergy = e;
    screenZ = a0/mom2;
    cosTetMaxNuc = std::max(cosThetaMax, 1.0 - 0.5*q2Limit/mom2);
    if(Z > 1.5) screenZ *= std::pow(Z,0.6666667)
		  *(1.13 + 3.76*invbeta2*Z*Z*chargeSquare*alpha2);
    else if(particle == theProton && cosTetMaxNuc < 0.0) 
      cosTetMaxNuc = 0.0;
    // A.V. Butkevich et al., NIM A 488 (2002) 282
    formfactA = mom2*constn*std::pow(A, 0.54);
  } 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
