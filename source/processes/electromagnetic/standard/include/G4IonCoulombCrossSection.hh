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
//	G4IonCoulombCrossSection.hh
//-------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:    G4IonCoulombCrossSection
//
// Author:      Cristina Consolandi
//
// Creation date: 05.10.2010 from G4eCoulombScatteringModel
//
// Class Description:
//      Computation of Screen-Coulomb Cross Section
//      for protons, alpha and heavy Ions
//
//
// Reference:
//      M.J. Boschini et al. "Nuclear and Non-Ionizing Energy-Loss
//      for Coulomb Scattered Particles from Low Energy up to Relativistic
//      Regime in Space Radiation Environment"
//      Accepted for publication in the Proceedings of  the  ICATPP Conference
//      on Cosmic Rays for Particle and Astroparticle Physics, Villa  Olmo, 7-8
//      October,  2010, to be published by World Scientific (Singapore).
//
//      Available for downloading at:
//      http://arxiv.org/abs/1011.4822
//
// -------------------------------------------------------------------

//
#ifndef G4IonCoulombCrossSection_h
#define G4IonCoulombCrossSection_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4NistManager.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4IonCoulombCrossSection
{
public:

  explicit G4IonCoulombCrossSection();

  ~G4IonCoulombCrossSection() = default;

  void Initialise(const G4ParticleDefinition*, G4double cosThetaLim);

  G4double NuclearCrossSection();

  G4double SampleCosineTheta();

  void SetupKinematic(G4double kinEnergy, G4double tmass);

  void SetupTarget(G4double Z, G4double kinEnergy, G4int heavycorr);

  inline void SetupParticle(const G4ParticleDefinition*);

  inline G4double GetMomentum2();

  G4IonCoulombCrossSection & operator=
  (const  G4IonCoulombCrossSection &right) = delete;
  G4IonCoulombCrossSection(const  G4IonCoulombCrossSection&) = delete;

private:

  void   SetScreenRSquare(G4int iz);

  const G4ParticleDefinition* theProton;  

  G4NistManager*  fNistManager;		
  G4Pow*          fG4pow;

  G4double                coeff;	  

  //cost - min - max 
  G4double                cosThetaMin;// def 1.0
  G4double                cosThetaMax;// def -1.0
  //SetupTarget
  G4double                cosTetMinNuc;// -->cosThetaMin
  G4double                cosTetMaxNuc;// -->cosThetaMax

  //cross section
  G4double                nucXSection;    	

  //energy 
  G4double                etag;	    

  // projectile........................
  const G4ParticleDefinition* particle;

  G4double                chargeSquare;  
  G4double                spin;	   
  G4double                mass;	   

  //lab of incedent particle 
  G4double                tkinLab;
  G4double                momLab2;
  G4double                invbetaLab2;

  //relative system with nucleus
  G4double                tkin;	   
  G4double                mom2;	   
  G4double                invbeta2;	   

  // target nucleus
  G4double                targetZ;    
  G4double                targetMass; 
  G4double                screenZ; 
  G4double                alpha2;
  G4double 	   	  ScreenRSquare;  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
void G4IonCoulombCrossSection::SetupParticle(const G4ParticleDefinition* p)
{
  particle = p;
  mass = particle->GetPDGMass();
  spin = particle->GetPDGSpin();
  if(0.0 != spin) { spin = 0.5; }
  G4double q = particle->GetPDGCharge()/CLHEP::eplus;
  chargeSquare = q*q;
  tkin = 0.0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4IonCoulombCrossSection::GetMomentum2()
{
  return mom2;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif


