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
// $Id: G4eCoulombScatteringModel.hh,v 1.8 2007-07-28 13:30:53 vnivanch Exp $
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

  G4double ScreeningParameter(G4double Z, G4double q2,
			      G4double mom2, G4double invbeta2);

  G4double NuclearSizeParameter(G4double A, G4double mom2);

  G4double CalculateECrossSectionPerAtom(const G4ParticleDefinition*, 
					 G4double kinEnergy, 
					 G4double Z, G4double ecut);

  G4double CalculateCrossSectionPerAtom(const G4ParticleDefinition*, 
					G4double kinEnergy, 
					G4double Z, G4double A);

private:

  // hide assignment operator
  G4eCoulombScatteringModel & operator=(const G4eCoulombScatteringModel &right);
  G4eCoulombScatteringModel(const  G4eCoulombScatteringModel&);

protected:

  G4ParticleChangeForGamma* fParticleChange;

  G4double                  coeff;
  G4double                  constn;
  G4double                  cosThetaMin;
  G4double                  cosThetaMax;
  G4double                  cosTetMaxNuc;
  G4double                  cosTetMaxElec;
  G4double                  q2Limit;
  G4double                  nucXS[100];
  G4double                  elXS[100];

private:

  const G4ParticleDefinition* theElectron;
  const G4ParticleDefinition* thePositron;
  const G4ParticleDefinition* theProton;

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
		G4double ecut, G4double)
{
  G4int iz = G4int(Z);
  G4bool b;
  if(theCrossSectionTable) {
    nucXS[iz] = std::exp((((*theCrossSectionTable)[index[G4int(Z)]]))
      ->GetValue(kinEnergy, b));
  } else nucXS[iz] = CalculateCrossSectionPerAtom(p, kinEnergy, Z, A);

  elXS[iz] = CalculateECrossSectionPerAtom(p, kinEnergy, Z, ecut);

  //  G4cout << "G4eCoulombScatteringModel:ComputeCSPerAtom e= " << kinEnergy 
  //	 << " Z= " << Z
  //       << "  CS= " << x << G4endl;
  return nucXS[iz] + elXS[iz];
}

inline G4double G4eCoulombScatteringModel::ScreeningParameter(
                G4double Z, 
                G4double q2, 
		G4double mom2, 
		G4double invbeta2)
{
  G4double R = a0/mom2;
  if(Z > 1.5) R *= std::pow(Z,0.6666667)*(1.13 + 3.76*invbeta2*Z*Z*q2*alpha2);
  return R;
}

inline G4double G4eCoulombScatteringModel::NuclearSizeParameter(
		G4double A, G4double mom2)
{
  // A.V. Butkevich et al., NIM A 488 (2002) 282
  return mom2*constn*std::pow(A, 0.54); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
