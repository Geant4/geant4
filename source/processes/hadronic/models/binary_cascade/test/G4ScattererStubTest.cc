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
#include "G4ScattererStub.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"


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

  G4ScattererStub scatterer;

  G4ParticleDefinition * proton = G4Proton::Proton();
  G4ParticleDefinition * neutron = G4Neutron::Neutron();
  G4double pmass = proton->GetPDGMass();
  G4double nmass = neutron->GetPDGMass();
  
  G4ThreeVector pos(0,0,0);

  for(G4int i = 0; i < nEvents; ++i)
  {
    G4double p = 1000*MeV*G4UniformRand();
    G4double theta = 2.0*G4UniformRand()-1.0;
    theta = acos(theta);
    G4double phi = G4UniformRand()*2*pi;
    G4ThreeVector direction(sin(theta)*cos(phi),
			    sin(theta)*sin(phi), cos(theta));
    G4LorentzVector mom1(p*direction, sqrt(p*p+pmass*pmass));
    p = 1000*MeV*G4UniformRand();
    theta = 2.0*G4UniformRand()-1.0;
    theta = acos(theta);
    phi = G4UniformRand()*2*pi;
    G4ThreeVector direction2(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    G4LorentzVector mom2(p*direction2, sqrt(p*p+nmass*nmass));

    G4KineticTrack kt1(proton, 0.0, pos, mom1);
    G4KineticTrack kt2(proton, 0.0, pos, mom2);

    G4KineticTrackVectorSTL * products = scatterer.Scatter(kt1, kt2);
    G4KineticTrackVectorSTL::iterator i = products->begin();
    G4LorentzVector fmom1 = (*i)->Get4Momentum();
    ++i;
    G4LorentzVector fmom2 = (*i)->Get4Momentum();
    G4LorentzVector diff = fmom1+fmom2-mom1-mom2;
    G4LorentzVector diff1 = fmom1-mom1;
    G4LorentzVector diff2 = fmom2-mom2;
    cout << diff1.x() << " " << diff1.y() << " " << diff1.z() << " " 
	 << diff1.t() << " " 
	 << diff2.x() << " " << diff2.y() << " " << diff2.z() << " " 
	 << diff2.t() << " " 
	 << diff.x() << " " << diff.y() << " " << diff.z() << " " 
	 << diff.t() << " " 
	 << fmom1.m() << " " << fmom2.m() << endl;
  }


  return 0;
}




