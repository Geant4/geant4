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
#include "G4Scatterer.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4IonConstructor.hh"


int main(int argc, char ** argv)
{
  G4int nEvents;
  if(argc == 2)
    nEvents = atoi(argv[1]);
  else
  {
    cout << "Input number of events: ";
    cin >> nEvents;
  }

 // constructing the particles
 
 G4LeptonConstructor aC1;
 G4BaryonConstructor aC2;
 G4MesonConstructor aC3;
 G4IonConstructor aC4;
 
 aC1.ConstructParticle();
 aC2.ConstructParticle();
 aC3.ConstructParticle();
 aC4.ConstructParticle();


  G4Scatterer scatterer;

  G4ParticleDefinition * proton = G4Proton::Proton();
  G4ParticleDefinition * neutron = G4Neutron::Neutron();
  G4double pmass = proton->GetPDGMass();
  G4double nmass = neutron->GetPDGMass();
  G4String pname=proton->GetParticleName();
  
  G4ThreeVector pos(0,0,0);

  for(G4int i = 0; i < nEvents; ++i)
  {
    G4double p = 2000*MeV *G4UniformRand();
    G4double theta = 2.0*G4UniformRand()-1.0;
    theta = acos(theta);
    G4double phi = G4UniformRand()*2*pi;
    G4ThreeVector direction(sin(theta)*cos(phi),
			    sin(theta)*sin(phi), cos(theta));
//    G4ThreeVector direction=G4ThreeVector(1.01745, -3.00003, -0.99999).unit();
    G4LorentzVector mom1(p*direction, sqrt(p*p+pmass*pmass));
    p = 200*MeV;		//*G4UniformRand();
    theta = 2.0*G4UniformRand()-1.0;
    theta = acos(theta);
    phi = G4UniformRand()*2*pi;
    G4ThreeVector direction2(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
//    G4ThreeVector direction2(0, -1, 0);
//  direction2 = -direction;

    G4LorentzVector mom2(p*direction2, sqrt(p*p+nmass*nmass));

    G4KineticTrack kt1(proton, 0.0, pos, mom1);
    G4KineticTrack kt2(neutron, 0.0, pos, mom2);

    G4KineticTrackVector * products = scatterer.Scatter(kt1, kt2);
    if ( products )
    {
	G4KineticTrackVector::iterator i = products->begin();
	G4LorentzVector fmom1 = (*i)->Get4Momentum();
	++i;
	G4LorentzVector fmom2 = (*i)->Get4Momentum();
	G4LorentzVector diff = fmom1+fmom2-mom1-mom2;
	G4LorentzVector diff1 = fmom1-mom1;
	G4LorentzVector diff2 = fmom2-mom2;
	G4cout << diff1<< " "
		<< diff2 << " "
		<< diff << " "
		<< fmom1.m() << " " << fmom2.m() << G4endl;
	G4cout << "in1/2, out 1/2 " << mom1 << " " << mom2 << " "
		<< fmom1 << " " << fmom2 << G4endl;
	G4cout << "costh particle 1 " << fmom1.vect().unit() * direction << G4endl;
    }
  }


  return 0;
}




