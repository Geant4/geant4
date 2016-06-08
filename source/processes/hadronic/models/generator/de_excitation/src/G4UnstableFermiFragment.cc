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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4UnstableFermiFragment.cc,v 1.7 2002/06/06 17:58:33 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4UnstableFermiFragment.hh"


G4UnstableFermiFragment::G4UnstableFermiFragment()
{
}

G4UnstableFermiFragment::G4UnstableFermiFragment(const G4UnstableFermiFragment &right)
{
  G4Exception("G4UnstableFermiFragment::copy_constructor meant to not be accessable");
}


G4UnstableFermiFragment::~G4UnstableFermiFragment()
{
}

  
const G4UnstableFermiFragment & G4UnstableFermiFragment::operator=(const G4UnstableFermiFragment &right)
{
  G4Exception("G4UnstableFermiFragment::operator= meant to not be accessable");
  return *this;
}


G4bool G4UnstableFermiFragment::operator==(const G4UnstableFermiFragment &right) const
{
  return false;
}

G4bool G4UnstableFermiFragment::operator!=(const G4UnstableFermiFragment &right) const
{
  return true;
}


G4std::deque<G4LorentzVector*> *
G4UnstableFermiFragment::FragmentsMomentum(G4double KinE, const G4int K, const G4double * Masses)
  // Calculates momentum for K fragments (Kopylov's method of sampling is used)
  // KinetEnergy is the available kinetic energy
{  
  G4std::deque<G4LorentzVector*>* MomentumList = 
    new G4std::deque<G4LorentzVector*>(K);


  G4double AvalaibleMass = 0; 
  for (G4int i=0; i<K; i++) AvalaibleMass += Masses[i];
  
  
  G4double PFragMagCM = 0.0;
  G4double Mass = AvalaibleMass+KinE;
  G4LorentzVector PFragCM(0.0,0.0,0.0,0.0);
  G4LorentzVector PFragLab(0.0,0.0,0.0,0.0);
  G4LorentzVector PRestCM(0.0,0.0,0.0,0.0);
  G4LorentzVector PRestLab(0.0,0.0,0.0,Mass);
  
  for (G4int l = 0; l < K-1; l++) {
    G4int LK = K - l;
  
    G4double FragMass = Masses[LK-1];

    AvalaibleMass -= FragMass;

    if (LK > 2) KinE *= RNKSI(LK-1);
    else KinE = 0.0;

    G4double RestMass = AvalaibleMass + KinE;


    PFragMagCM = sqrt(
		      abs((Mass*Mass - (FragMass + RestMass)*(FragMass + RestMass))*
			  (Mass*Mass - (FragMass - RestMass)*(FragMass - RestMass)))
		      )/(2.0*Mass);


    // Create a unit vector with a random direction isotropically distributed
    G4ParticleMomentum RandVector(IsotropicVector(PFragMagCM));

    PFragCM.setVect(RandVector);
    PFragCM.setE(sqrt(RandVector.mag2()+FragMass*FragMass));

    PRestCM.setVect(-RandVector);
    PRestCM.setE(sqrt(RandVector.mag2()+RestMass*RestMass));

    

    G4ThreeVector BoostV = PRestLab.boostVector();

    PFragLab = PFragCM;
    PFragLab.boost(BoostV);

    PRestLab = PRestCM;
    PRestLab.boost(BoostV);


    MomentumList->push_front(new G4LorentzVector(PFragLab));

    Mass = RestMass;
  }

  MomentumList->push_front(new G4LorentzVector(PRestLab));
  return MomentumList;
}


G4double G4UnstableFermiFragment::RNKSI(const G4int K)
{
  G4double csim = (3.0*K-5.0)/(3.0*K-4.0);
  G4double pex = 1.5*K-2.5;
  G4double fcsim = sqrt(1.0-csim)*pow(csim,pex);

  G4double csi = 0.0;
  G4double fcsi= 0.0;
  G4double rf = 0.0;
  do {
    csi = G4UniformRand();
    fcsi = sqrt(1.0-csi)*pow(csi,pex);
    rf = fcsim*G4UniformRand();
  } while (rf > fcsi);
  return csi;
}
    



G4ParticleMomentum G4UnstableFermiFragment::IsotropicVector(const G4double Magnitude)
  // Samples a isotropic random vectorwith a magnitud given by Magnitude.
  // By default Magnitude = 1.0
{
  G4double CosTheta = 1.0 - 2.0*G4UniformRand();
  G4double SinTheta = sqrt(1.0 - CosTheta*CosTheta);
  G4double Phi = twopi*G4UniformRand();
  G4ParticleMomentum Vector(Magnitude*cos(Phi)*SinTheta,
			    Magnitude*sin(Phi)*SinTheta,
			    Magnitude*CosTheta);
  return Vector;
}
