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
// $Id: G4WentzelOKandVIxSection.cc 98737 2016-08-09 12:51:38Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4WentzelOKandVIxSection
//
// Author:      V.Ivanchenko 
//
// Creation date: 09.04.2008 from G4MuMscModel
//
// Modifications:
//
// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4WentzelOKandVIxSection.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4EmParameters.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelOKandVIxSection::ScreenRSquareElec[] = {0.0};
G4double G4WentzelOKandVIxSection::ScreenRSquare[]     = {0.0};
G4double G4WentzelOKandVIxSection::FormFactor[]        = {0.0};

using namespace std;

G4WentzelOKandVIxSection::G4WentzelOKandVIxSection(G4bool combined) :
  temp(0.,0.,0.),
  numlimit(0.1),
  nwarnings(0),
  nwarnlimit(50),
  isCombined(combined),
  cosThetaMax(-1.0),
  alpha2(fine_structure_const*fine_structure_const)
{
  fNistManager = G4NistManager::Instance();
  fG4pow = G4Pow::GetInstance();
  theElectron = G4Electron::Electron();
  thePositron = G4Positron::Positron();
  theProton   = G4Proton::Proton();
  lowEnergyLimit = 1.0*eV;
  G4double p0 = electron_mass_c2*classic_electr_radius;
  coeff = twopi*p0*p0;
  particle = nullptr;

  fNucFormfactor = fExponentialNF;

  // Thomas-Fermi screening radii
  // Formfactors from A.V. Butkevich et al., NIM A 488 (2002) 282

  if(0.0 == ScreenRSquare[0]) {
    G4double a0 = electron_mass_c2/0.88534; 
    G4double constn = 6.937e-6/(MeV*MeV);

    ScreenRSquare[0] = alpha2*a0*a0;
    ScreenRSquareElec[0] = ScreenRSquare[0];
    for(G4int j=1; j<100; ++j) {
      G4double x = a0*fG4pow->Z13(j);
      if(1 == j) { ScreenRSquare[j] = 0.5*alpha2*a0*a0; }
      else {
        ScreenRSquare[j] = 0.5*(1 + G4Exp(-j*j*0.001))*alpha2*x*x;
        ScreenRSquareElec[j] = 0.5*alpha2*x*x;
      }
      x = fNistManager->GetA27(j);
      FormFactor[j] = constn*x*x;
    } 
  }
  currentMaterial = 0;
  factB = factD = formfactA = screenZ = 0.0;
  cosTetMaxElec = cosTetMaxNuc = invbeta2 = kinFactor = gam0pcmp = pcmp2 = 1.0;

  factB1= 0.5*CLHEP::pi*fine_structure_const;

  tkin = mom2 = momCM2 = factorA2 = mass = spin = chargeSquare = charge3 = 0.0;
  ecut = etag = DBL_MAX;
  targetZ = 0;
  targetMass = proton_mass_c2;
  particle = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4WentzelOKandVIxSection::~G4WentzelOKandVIxSection()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4WentzelOKandVIxSection::Initialise(const G4ParticleDefinition* p, 
                                          G4double cosThetaLim)
{
  SetupParticle(p);
  tkin = mom2 = momCM2 = 0.0;
  ecut = etag = DBL_MAX;
  targetZ = 0;

  // cosThetaMax is below 1.0 only when MSC is combined with SS
  if(isCombined) { cosThetaMax = cosThetaLim; }
  
  G4double a = G4EmParameters::Instance()->FactorForAngleLimit()
    *CLHEP::hbarc/CLHEP::fermi;
  factorA2 = 0.5*a*a;
  currentMaterial = nullptr;

  fNucFormfactor = G4EmParameters::Instance()->NuclearFormfactorType();
  
  //G4cout << "G4WentzelOKandVIxSection::Initialise  mass= " << mass
  //         << "  " << p->GetParticleName() 
  //     << "  cosThetaMax= " << cosThetaMax << G4endl; 
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4WentzelOKandVIxSection::SetupParticle(const G4ParticleDefinition* p)
{
  particle = p;
  mass = particle->GetPDGMass();
  spin = particle->GetPDGSpin();
  if(0.0 != spin) { spin = 0.5; }
  G4double q = std::fabs(particle->GetPDGCharge()/eplus);
  chargeSquare = q*q;
  charge3 = chargeSquare*q;
  tkin = 0.0;
  currentMaterial = 0;
  targetZ = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
G4double
G4WentzelOKandVIxSection::SetupTarget(G4int Z, G4double cut)
{
  G4double cosTetMaxNuc2 = cosTetMaxNuc;
  if(Z != targetZ || tkin != etag) {
    etag    = tkin; 
    targetZ = Z;
    if(targetZ > 99) { targetZ = 99; }
    G4double massT = proton_mass_c2;
    if(targetZ > 1) { 
      massT = fNistManager->GetAtomicMassAmu(targetZ)*CLHEP::amu_c2;
    }
    SetTargetMass(massT);

    kinFactor = coeff*Z*chargeSquare*invbeta2/mom2;

    if(1 == Z) {
      screenZ = ScreenRSquare[targetZ]/mom2;
    } else if(mass > MeV) {
      screenZ = std::min(Z*1.13,1.13 +3.76*Z*Z*invbeta2*alpha2*chargeSquare)*
        ScreenRSquare[targetZ]/mom2;
    } else {
      G4double tau = tkin/mass;
      screenZ = std::min(Z*1.13,(1.13 +3.76*Z*Z
          *invbeta2*alpha2*std::sqrt(tau/(tau + fG4pow->Z23(targetZ)))))*
        ScreenRSquareElec[targetZ]/mom2;
    }
    if(targetZ == 1 && cosTetMaxNuc2 < 0.0 && particle == theProton) {
      cosTetMaxNuc2 = 0.0;
    }
    formfactA = FormFactor[targetZ]*mom2;

    cosTetMaxElec = 1.0;
    ComputeMaxElectronScattering(cut); 
  }
  return cosTetMaxNuc2;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4WentzelOKandVIxSection::ComputeTransportCrossSectionPerAtom(G4double cosTMax)
{
  G4double xSection = 0.0;
  if(cosTMax >= 1.0) { return xSection; }
 
  G4double x = 0; 
  G4double y = 0;
  G4double x1= 0;
  G4double x2= 0;
  G4double xlog = 0.0;

  G4double costm = std::max(cosTMax,cosTetMaxElec); 
  G4double fb = screenZ*factB;

  // scattering off electrons
  if(costm < 1.0) {
    x = (1.0 - costm)/screenZ;
    if(x < numlimit) { 
      x2 = 0.5*x*x;
      y  = x2*(1.0 - 1.3333333*x + 3*x2);
      if(0.0 < factB) { y -= fb*x2*x*(0.6666667 - x); }
    } else { 
      x1= x/(1 + x);
      xlog = G4Log(1.0 + x);  
      y = xlog - x1; 
      if(0.0 < factB) { y -= fb*(x + x1 - 2*xlog); }
    }

    if(y < 0.0) {
      ++nwarnings;
      if(nwarnings < nwarnlimit) {
        G4cout << "G4WentzelOKandVIxSection::ComputeTransportCrossSectionPerAtom"
               << " scattering on e- <0"
               << G4endl;
        G4cout << "y= " << y 
               << " e(MeV)= " << tkin << " p(MeV/c)= " << sqrt(mom2) 
               << " Z= " << targetZ << "  " 
               << particle->GetParticleName() << G4endl;
        G4cout << " 1-costm= " << 1.0-costm << " screenZ= " << screenZ 
               << " x= " << x << G4endl;
      }
      y = 0.0;
    }
    xSection = y;
  }
  /* 
       G4cout << "G4WentzelVI:XS per A " << " Z= " << targetZ 
       << " e(MeV)= " << tkin/MeV << " XSel= " << xSection
         << " cut(MeV)= " << ecut/MeV  
           << " zmaxE= " << (1.0 - cosTetMaxElec)/screenZ 
         << " zmaxN= " << (1.0 - cosThetaMax)/screenZ 
         << " 1-costm= " << 1.0 - cosThetaMax << G4endl;
  */
  // scattering off nucleus
  if(cosTMax < 1.0) {
    x = (1.0 - cosTMax)/screenZ;
    if(x < numlimit) { 
      x2 = 0.5*x*x;
      y  = x2*(1.0 - 1.3333333*x + 3*x2); 
      if(0.0 < factB) { y -= fb*x2*x*(0.6666667 - x); }
    } else { 
      x1= x/(1 + x);
      xlog = G4Log(1.0 + x);  
      y = xlog - x1; 
      if(0.0 < factB) { y -= fb*(x + x1 - 2*xlog); }
    }

    if(y < 0.0) {
      ++nwarnings;
      if(nwarnings < nwarnlimit) {
        G4cout << "G4WentzelOKandVIxSection::ComputeTransportCrossSectionPerAtom"
               << " scattering on nucleus <0"
               << G4endl;
        G4cout << "y= " << y 
               << " e(MeV)= " << tkin << " Z= " << targetZ << "  " 
               << particle->GetParticleName() << G4endl;
        G4cout << " formfactA= " << formfactA << " screenZ= " << screenZ 
               << " x= " << " x1= " << x1 <<G4endl;
      }
      y = 0.0;
    }
    xSection += y*targetZ; 
  }
  xSection *= kinFactor;
 
  /*
  G4cout << "Z= " << targetZ << " XStot= " << xSection/barn 
         << " screenZ= " << screenZ << " formF= " << formfactA 
         << " for " << particle->GetParticleName() 
           << " m= " << mass << " 1/v= " << sqrt(invbeta2) << " p= " << sqrt(mom2)
         << " x= " << x 
         << G4endl;
  */
  return xSection; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector&
G4WentzelOKandVIxSection::SampleSingleScattering(G4double cosTMin,
                                                 G4double cosTMax,
                                                 G4double elecRatio)
{
  temp.set(0.0,0.0,1.0);
  CLHEP::HepRandomEngine* rndmEngineMod = G4Random::getTheEngine();
 
  G4double formf = formfactA;
  G4double cost1 = cosTMin;
  G4double cost2 = cosTMax;
  if(elecRatio > 0.0) {
    if(rndmEngineMod->flat() <= elecRatio) {
      formf = 0.0;
      cost1 = std::max(cost1,cosTetMaxElec);
      cost2 = std::max(cost2,cosTetMaxElec);
    }
  }
  if(cost1 > cost2) {
    G4double w1 = 1. - cost1 + screenZ;
    G4double w2 = 1. - cost2 + screenZ;
    G4double z1 = w1*w2/(w1 + rndmEngineMod->flat()*(w2 - w1)) - screenZ;

    G4double fm = 1.0;
    if(fNucFormfactor == fExponentialNF) {
      fm += formf*z1;
      fm = 1.0/(fm*fm);
    } else if(fNucFormfactor == fGaussianNF) {
      fm = G4Exp(-2*formf*z1);
    } else if(fNucFormfactor == fFlatNF) {
      static const G4double ccoef = 0.00508/MeV;
      G4double x = std::sqrt(2.*mom2*z1)*ccoef*2.;
      fm = FlatFormfactor(x);
      fm *= FlatFormfactor(x*0.6
	    *fG4pow->A13(fNistManager->GetAtomicMassAmu(targetZ) ) );     
    }

    G4double grej = (1. - z1*factB + factB1*targetZ*sqrt(z1*factB)*(2. - z1))
      *fm*fm/(1.0 + z1*factD);

    if(rndmEngineMod->flat() <= grej) {
      // exclude "false" scattering due to formfactor and spin effect
      G4double cost = 1.0 - z1;
      if(cost > 1.0)       { cost = 1.0; }
      else if(cost < -1.0) { cost =-1.0; }
      G4double sint = sqrt((1.0 - cost)*(1.0 + cost));
      //G4cout << "sint= " << sint << G4endl;
      G4double phi  = twopi*rndmEngineMod->flat();
      temp.set(sint*cos(phi),sint*sin(phi),cost);
    }
  }
  return temp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4WentzelOKandVIxSection::ComputeMaxElectronScattering(G4double cutEnergy)
{
  if(mass > MeV) {
    G4double ratio = electron_mass_c2/mass;
    G4double tau = tkin/mass;
    G4double tmax = 2.0*electron_mass_c2*tau*(tau + 2.)/
      (1.0 + 2.0*ratio*(tau + 1.0) + ratio*ratio);
    //tmax = std::min(tmax, targetZ*targetZ*10*eV); 
    cosTetMaxElec = 1.0 - std::min(cutEnergy, tmax)*electron_mass_c2/mom2;
  } else {

    G4double tmax = tkin;
    if(particle == theElectron) { tmax *= 0.5; }
    //tmax = std::min(tmax, targetZ*targetZ*10*eV); 
    G4double t = std::min(cutEnergy, tmax);
    G4double mom21 = t*(t + 2.0*electron_mass_c2);
    G4double t1 = tkin - t;
    //G4cout <<"tkin=" <<tkin<<" tmax= "<<tmax<<" t= " 
    //<<t<< " t1= "<<t1<<" cut= "<<ecut<<G4endl;
    if(t1 > 0.0) {
      G4double mom22 = t1*(t1 + 2.0*mass);
      G4double ctm = (mom2 + mom22 - mom21)*0.5/sqrt(mom2*mom22);
      if(ctm <  1.0) { cosTetMaxElec = ctm; }
      if(particle == theElectron && cosTetMaxElec < 0.0) { 
        cosTetMaxElec = 0.0; 
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4WentzelOKandVIxSection::ComputeSecondTransportMoment(G4double /*CosThetaMax*/)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
