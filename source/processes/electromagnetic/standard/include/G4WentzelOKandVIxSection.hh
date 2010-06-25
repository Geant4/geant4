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
// $Id: G4WentzelOKandVIxSection.hh,v 1.6 2010-06-25 09:41:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4WentzelOKandVIxSection
//
// Authors:       V.Ivanchenko and O.Kadri 
//
// Creation date: 21.05.2010 
//
// Modifications:
//
//
// Class Description:
//
// Implementation of the computation of total and transport cross sections,
// sample scattering angle for the single scattering case.
// to be used by single and multiple scattering models. References:
// 1) G.Wentzel, Z. Phys. 40 (1927) 590.
// 2) J.M. Fernandez-Varea et al., NIM B73 (1993) 447.
//
// -------------------------------------------------------------------
//

#ifndef G4WentzelOKandVIxSection_h
#define G4WentzelOKandVIxSection_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4NistManager.hh"
#include "G4ThreeVector.hh"
#include "G4Pow.hh"

class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4WentzelOKandVIxSection 
{

public:

  G4WentzelOKandVIxSection();

  virtual ~G4WentzelOKandVIxSection();

  void Initialise(const G4ParticleDefinition*, G4double CosThetaLim);

  void SetupParticle(const G4ParticleDefinition*);

  // return cos(ThetaMax) for msc and cos(thetaMin) for single scattering
  // cut = DBL_MAX means no scattering off electrons 
  G4double SetupTarget(G4int Z, G4double cut = DBL_MAX);

  G4double ComputeTransportCrossSectionPerAtom(G4double CosThetaMax);
 
  G4ThreeVector SampleSingleScattering(G4double CosThetaMin,
				       G4double CosThetaMax,
				       G4double elecRatio = 0.0);

  inline G4double ComputeNuclearCrossSection(G4double CosThetaMin,
					     G4double CosThetaMax);
 
  inline G4double ComputeElectronCrossSection(G4double CosThetaMin,
					      G4double CosThetaMax);
 
  inline G4double SetupKinematic(G4double kinEnergy, const G4Material* mat);
  
  inline void SetTargetMass(G4double value);

  inline void SetRelativisticMass(G4double value);

  inline G4double GetMomentumSquare() const;

  inline G4double GetCosThetaNuc() const;

  inline G4double GetCosThetaElec() const;

private:

  void ComputeMaxElectronScattering(G4double cut);

  //  hide assignment operator
  G4WentzelOKandVIxSection & operator=(const  G4WentzelOKandVIxSection &right);
  G4WentzelOKandVIxSection(const  G4WentzelOKandVIxSection&);

  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theElectron;
  const G4ParticleDefinition* thePositron;
  const G4Material* currentMaterial;

  G4NistManager*  fNistManager;
  G4Pow*          fG4pow;

  G4double numlimit;

  G4double elecXSRatio;

  // integer parameters
  G4int    nwarnings;
  G4int    nwarnlimit;

  // single scattering parameters
  G4double coeff;
  G4double cosTetMaxElec;
  G4double cosTetMaxNuc;
  G4double cosThetaMax;
  G4double alpha2;

  // projectile
  const G4ParticleDefinition* particle;

  G4double chargeSquare;
  G4double charge3;
  G4double spin;
  G4double mass;
  G4double tkin;
  G4double mom2;
  G4double invbeta2;
  G4double kinFactor;
  G4double etag;
  G4double ecut;
  G4double lowEnergyLimit;

  // target
  G4int    targetZ;
  G4double targetMass;
  G4double screenZ;
  G4double formfactA;
  G4double factorA2;
  G4double factB;
  G4double factD;

  static G4double ScreenRSquare[100];
  static G4double FormFactor[100];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4WentzelOKandVIxSection::SetupKinematic(G4double ekin, const G4Material* mat)
{
  if(ekin != tkin || mat != currentMaterial) { 
    currentMaterial = mat;
    tkin  = ekin;
    mom2  = tkin*(tkin + 2.0*mass);
    invbeta2 = 1.0 +  mass*mass/mom2;
    factB = spin/invbeta2;
    cosTetMaxNuc = cosThetaMax;
    if(std::fabs(cosThetaMax) < 1.0) {
      cosTetMaxNuc = 
	std::max(cosThetaMax,1.-factorA2*mat->GetIonisation()->GetInvA23()/mom2);
    }
  } 
  return cosTetMaxNuc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4WentzelOKandVIxSection::SetTargetMass(G4double value)
{
  targetMass = value;
  factD = std::sqrt(mom2)/value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4WentzelOKandVIxSection::SetRelativisticMass(G4double value)
{
  mass = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4WentzelOKandVIxSection::GetMomentumSquare() const
{
  return mom2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4WentzelOKandVIxSection::GetCosThetaNuc() const
{
  return cosTetMaxNuc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4WentzelOKandVIxSection::GetCosThetaElec() const
{
  return cosTetMaxElec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
G4WentzelOKandVIxSection::ComputeNuclearCrossSection(G4double cosTMin,
						     G4double cosTMax)
{
  G4double xsec = 0.0;
  if(cosTMax < cosTMin) {
    xsec = targetZ*kinFactor*(cosTMin - cosTMax)/
      ((1.0 - cosTMin + screenZ)*(1.0 - cosTMax + screenZ));
  }
  return xsec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
G4WentzelOKandVIxSection::ComputeElectronCrossSection(G4double cosTMin,
						      G4double cosTMax)
{
  G4double xsec = 0.0;
  G4double cost1 = std::max(cosTMin,cosTetMaxElec);
  G4double cost2 = std::max(cosTMax,cosTetMaxElec);
  if(cost1 > cost2) {
    xsec = kinFactor*(cost1 - cost2)/((1.0 - cost1 + screenZ)*(1.0 - cost2 + screenZ));
  }
  return xsec;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

