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
#include "globals.hh"
#include "G4Absorber.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4KineticTrackVectorSTL.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"
#include "G4PionPlus.hh"
#include "G4LeptonConstructor.hh"


G4V3DNucleus * the3DNucleus;
G4KineticTrackVectorSTL theTargetList;

void BuildTargetList();
G4ThreeVector GetSpherePoint(G4double r);

int main(int argc, char ** argv)
{
  G4LeptonConstructor leptons;
  leptons.ConstructParticle();

  G4int A, Z;
  G4int nEvents;
  if(argc != 4)
  {
    cout << "Input A and Z: ";
    cin >> A >> Z;
    cout << "Input number of events: ";
    cin >> nEvents;
  }
  else
  {
    A = atoi(argv[1]);
    Z = atoi(argv[2]);
    nEvents = atoi(argv[3]);
  }

// create the nucleus
  the3DNucleus = new G4Fancy3DNucleus;
  the3DNucleus->Init(A, Z);
// create theTargetList
  BuildTargetList();
// create the pion
    G4ParticleDefinition * pion = G4PionPlus::PionPlus();
    G4ThreeVector pos = GetSpherePoint(the3DNucleus->GetOuterRadius());
    G4double p = 100*MeV;
    G4double mass = pion->GetPDGMass();
    G4LorentzVector mom(0, 0, p, sqrt(p*p+mass*mass));
    G4KineticTrack kt(pion, 0., pos, mom);
// create the absorber
    G4double theCutOnP = 150*MeV;
    G4Absorber absorber(theCutOnP);

    absorber.FindAbsorbers(kt, theTargetList);
  for(G4int i = 0; i < nEvents; ++i)
  {
    absorber.FindProducts(kt);
  }
  
  return 0;
}


void BuildTargetList()
{
  if(!the3DNucleus->StartLoop())
  {
    G4cerr << "G4HadronKineticModel::BuildTargetList(): StartLoop() error!"
	   << G4endl;
    return;
  }
  G4Nucleon * nucleon;
  G4ParticleDefinition * definition;
  G4ThreeVector pos;
  G4LorentzVector mom;
  while((nucleon = the3DNucleus->GetNextNucleon()) != NULL)
  {
    definition = nucleon->GetDefinition();
    pos = nucleon->GetPosition();
    mom = nucleon->GetMomentum();
    G4KineticTrack * kt = new G4KineticTrack(definition, 0., pos, mom);
    theTargetList.push_back(kt);
  }
}


G4ThreeVector GetSpherePoint(G4double r)
{
// Get the entry point of the projectile, distribuited uniformly
// on the projection of the surface in the plane ortogonal to the direction
// of the motion.
  G4double b = r*G4UniformRand();  // impact parameter
  G4double phi = G4UniformRand()*2*pi;
  G4double x = b*cos(phi);
  G4double y = b*sin(phi);
  G4double z = -sqrt(r*r-b*b);
  z *= 1.001; // Get position a little bit out of the sphere...
  G4ThreeVector point;
  point.setX(x);
  point.setY(y);
  point.setZ(z);
  return point;
}
