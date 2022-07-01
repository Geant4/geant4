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
//      G4IonCoulombCrossSection.cc
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
//	Computation of Screen-Coulomb Cross Section 
//	for protons, alpha and heavy Ions
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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4IonCoulombCrossSection.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4Exp.hh"
#include "G4Log.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

const G4double a0 = CLHEP::electron_mass_c2/0.88534;

G4IonCoulombCrossSection::G4IonCoulombCrossSection():
   cosThetaMin(1.0),
   cosThetaMax(-1.0),
   alpha2(fine_structure_const*fine_structure_const)
{
  fNistManager = G4NistManager::Instance();
  fG4pow = G4Pow::GetInstance();
  theProton   = G4Proton::Proton();
  particle = nullptr;

  G4double p0 = electron_mass_c2*classic_electr_radius;
  coeff  = twopi*p0*p0;

  cosTetMinNuc=0;
  cosTetMaxNuc=0;
  nucXSection =0;

  chargeSquare = spin = mass = 0.0;
  tkinLab = momLab2 = invbetaLab2 = tkin = mom2 = invbeta2 = 0.0;

  targetZ = targetMass = screenZ = ScreenRSquare = etag = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4IonCoulombCrossSection::Initialise(const G4ParticleDefinition* p,
                                          G4double CosThetaLim)
{
  SetupParticle(p);
  nucXSection = tkin = targetZ = mom2 = 0.0;
  etag = DBL_MAX;
  particle = p;		
  cosThetaMin = CosThetaLim; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4IonCoulombCrossSection::SetupKinematic(G4double ekin, G4double tmass)
{
  if(ekin != tkinLab || tmass != targetMass) {

    // lab
    tkinLab = ekin;
    momLab2 = tkinLab*(tkinLab + 2.0*mass);
    invbetaLab2 = 1.0 +  mass*mass/momLab2;

    G4double etot = tkinLab + mass;
    G4double ptot = sqrt(momLab2);
    G4double m12  = mass*mass;
    // relativistic reduced mass from publucation
    // A.P. Martynenko, R.N. Faustov, Teoret. mat. Fiz. 64 (1985) 179
        
    //incident particle & target nucleus
    targetMass = tmass;
    G4double Ecm=sqrt(m12 + targetMass*targetMass + 2.0*etot*targetMass);
    G4double mu_rel=mass*targetMass/Ecm;
    G4double momCM= ptot*targetMass/Ecm;
    // relative system
    mom2 = momCM*momCM;
    invbeta2 = 1.0 +  mu_rel*mu_rel/mom2;
    tkin = momCM*sqrt(invbeta2) - mu_rel;//Ekin of mu_rel

    cosTetMinNuc = cosThetaMin;
    cosTetMaxNuc = cosThetaMax;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4IonCoulombCrossSection::SetupTarget(G4double Z, G4double e, 
					   G4int)
{
  if(Z != targetZ || e != etag) {
    etag    = e;
    targetZ = Z;
    G4int iz= G4lrint(Z);

    SetScreenRSquare(iz);
    screenZ = 0;
    screenZ = ScreenRSquare/mom2;
    //heavycorr = 0;
    //	G4cout<< "heavycorr "<<heavycorr<<G4endl;

    G4double corr=5.*twopi*Z*std::sqrt(chargeSquare*alpha2);
    corr=G4Exp(G4Log(corr)*0.04);
    screenZ *=0.5*(1.13 + corr*3.76*Z*Z*chargeSquare*invbeta2*alpha2);
    // G4cout<<" heavycorr Z e corr....2As "<< heavycorr << "\t"
    //  <<Z <<"\t"<<e/MeV <<"\t"<<screenZ<<G4endl;
      
    if(1 == iz && particle == theProton && cosTetMaxNuc < 0.0) {
      cosTetMaxNuc = 0.0;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4IonCoulombCrossSection::SetScreenRSquare(G4int iz)
{
  //for proton Thomas-Fermi screening length    
  G4int Z1 = G4lrint(std::sqrt(chargeSquare));
  G4double Z113  = fG4pow->Z13(iz);
  G4double Z1023 = fG4pow->powZ(Z1,0.23);
  G4double Z2023 = fG4pow->powZ(iz,0.23);
  G4double x=a0*(Z1023+Z2023);  
              
  // Universal screening length
  if(particle == theProton){
     x = a0*Z113;
  } 

  ScreenRSquare = alpha2*x*x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonCoulombCrossSection::NuclearCrossSection()
{
  // This method needs initialisation before be called
  // scattering with target nucleus
  G4double fac = coeff*targetZ*(targetZ)*chargeSquare*invbeta2/mom2;

  nucXSection  = 0.0;

  G4double x  = 1.0 - cosTetMinNuc;
  G4double x1 = x + screenZ;

  // scattering with nucleus
  if(cosTetMaxNuc < cosTetMinNuc) {
    nucXSection = fac*(cosTetMinNuc - cosTetMaxNuc)/
      (x1*(1.0 - cosTetMaxNuc + screenZ));
  }

  return nucXSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4IonCoulombCrossSection::SampleCosineTheta()
{  	
  G4double z1 = 0.0;
  if(cosTetMaxNuc < cosTetMinNuc) {

    G4double x1 = 1. - cosTetMinNuc + screenZ;
    G4double x2 = 1. - cosTetMaxNuc + screenZ;
    G4double dx = cosTetMinNuc - cosTetMaxNuc;
    z1 = x1*x2/(x1 + G4UniformRand()*dx) - screenZ;
  }
  return z1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




