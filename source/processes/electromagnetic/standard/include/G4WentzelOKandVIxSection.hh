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
#include "G4NuclearFormfactorType.hh"
#include "G4ThreeVector.hh"
#include "G4Pow.hh"

class G4ParticleDefinition;
class G4ScreeningMottCrossSection;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4WentzelOKandVIxSection 
{

public:

  explicit G4WentzelOKandVIxSection(G4bool comb=true);

  virtual ~G4WentzelOKandVIxSection();

  void Initialise(const G4ParticleDefinition*, G4double CosThetaLim);

  void SetupParticle(const G4ParticleDefinition*);

  // return cos(ThetaMax) for msc and cos(thetaMin) for single scattering
  // cut = DBL_MAX means no scattering off electrons 
  virtual G4double SetupKinematic(G4double kinEnergy, const G4Material* mat);
  G4double SetupTarget(G4int Z, G4double cut);

  G4double ComputeTransportCrossSectionPerAtom(G4double CosThetaMax);
 
  G4ThreeVector& SampleSingleScattering(G4double CosThetaMin,
					G4double CosThetaMax,
					G4double elecRatio);

  G4double ComputeSecondTransportMoment(G4double CosThetaMax);

  inline G4double ComputeNuclearCrossSection(G4double CosThetaMin,
					     G4double CosThetaMax);
 
  inline G4double ComputeElectronCrossSection(G4double CosThetaMin,
					      G4double CosThetaMax);
   
  inline void SetTargetMass(G4double value);

  inline G4double GetMomentumSquare() const;

  inline G4double GetCosThetaNuc() const;

  inline G4double GetCosThetaElec() const;

  //  hide assignment operator
  G4WentzelOKandVIxSection & operator=
  (const G4WentzelOKandVIxSection &right) = delete;
  G4WentzelOKandVIxSection(const  G4WentzelOKandVIxSection&) = delete;

protected:

  void ComputeMaxElectronScattering(G4double cut);

  void InitialiseA();

  inline G4double FlatFormfactor(G4double x);

  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theElectron;
  const G4ParticleDefinition* thePositron;
  const G4ParticleDefinition* particle = nullptr;
  const G4Material* currentMaterial = nullptr;

  G4NistManager* fNistManager;
  G4Pow*         fG4pow;

  G4ScreeningMottCrossSection* fMottXSection = nullptr;

  G4ThreeVector temp;

  // single scattering parameters
  G4double coeff;
  G4double cosTetMaxElec = 1.0;
  G4double cosTetMaxNuc = 1.0;
  G4double cosThetaMax = -1.0;

  G4double chargeSquare = 0.0;
  G4double charge3 = 0.0;
  G4double spin = 0.0;
  G4double mass = 0.0;
  G4double tkin = 0.0;
  G4double mom2 = 0.0;
  G4double momCM2 = 0.0;
  G4double invbeta2 = 1.0;
  G4double kinFactor = 1.0;
  G4double etag = DBL_MAX;
  G4double ecut = DBL_MAX;

  // target
  G4double targetMass;
  G4double screenZ = 0.0;
  G4double formfactA = 0.0;
  G4double factorA2 = 0.0;
  G4double factB = 0.0;
  G4double factD = 0.0;
  G4double fMottFactor = 1.0;
  G4double gam0pcmp = 1.0;
  G4double pcmp2 = 1.0;

  // integer parameters
  G4int targetZ = 0;
  G4int nwarnings = 0;

  G4NuclearFormfactorType fNucFormfactor = fExponentialNF;

  G4bool isCombined;

  static G4double ScreenRSquareElec[100];
  static G4double ScreenRSquare[100];
  static G4double FormFactor[100];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4WentzelOKandVIxSection::SetTargetMass(G4double value)
{
  targetMass = value;
  factD = std::sqrt(mom2)/value;
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
  return targetZ*kinFactor*fMottFactor*(cosTMin - cosTMax)/
    ((1.0 - cosTMin + screenZ)*(1.0 - cosTMax + screenZ));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double 
G4WentzelOKandVIxSection::ComputeElectronCrossSection(G4double cosTMin,
						      G4double cosTMax)
{
  G4double cost1 = std::max(cosTMin,cosTetMaxElec);
  G4double cost2 = std::max(cosTMax,cosTetMaxElec);
  return (cost1 <= cost2) ? 0.0 : kinFactor*fMottFactor*(cost1 - cost2)/
    ((1.0 - cost1 + screenZ)*(1.0 - cost2 + screenZ));
}

inline G4double G4WentzelOKandVIxSection::FlatFormfactor(G4double x)
{
  return 3.0*(std::sin(x) - x*std::cos(x))/(x*x*x);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

