// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4CollisionTest.cc
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it), 
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"

#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"
#include "g4std/iostream"
#include <assert.h>

#include "Randomize.hh"

#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"

#include "G4AntiProton.hh"
#include "G4Proton.hh"
#include "G4AntiNeutron.hh"
#include "G4Neutron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4KaonMinus.hh"
#include "G4KaonPlus.hh"
#include "G4Gamma.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4Lambda.hh"

#include "G4VShortLivedParticle.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedTable.hh"

#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4KineticTrackVectorSTL.hh"

#include "G4ParticleDefinition.hh"

#include "G4VCollision.hh"
#include "G4CollisionNN.hh"
#include "G4CollisionMesonBaryon.hh"

int main()
{

  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();

  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  G4ParticleDefinition* antiProton = G4AntiProton::AntiProtonDefinition();
  G4ParticleDefinition* neutron = G4Neutron::NeutronDefinition();   
  G4ParticleDefinition* antiNeutron = G4AntiNeutron::AntiNeutronDefinition();   

  G4ParticleDefinition* pionPlus = G4PionPlus::PionPlusDefinition();
  G4ParticleDefinition* pionMinus = G4PionMinus::PionMinusDefinition();
  G4ParticleDefinition* pionZero = G4PionZero::PionZeroDefinition();

  G4ParticleDefinition* kaonPlus = G4KaonPlus::KaonPlusDefinition();
  G4ParticleDefinition* kaonMinus = G4KaonMinus::KaonMinusDefinition();
   
  G4ParticleDefinition* lambda = G4Lambda::LambdaDefinition();

  G4ParticleDefinition* theMuonPlus = G4MuonPlus::MuonPlusDefinition();
  G4ParticleDefinition* theMuonMinus = G4MuonMinus::MuonMinusDefinition();

  G4ParticleDefinition* theNeutrinoMu = G4NeutrinoMu::NeutrinoMuDefinition();
  G4ParticleDefinition* theAntiNeutrinoMu = G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();

  // Select track 1
  

  G4ParticleDefinition* primary1 = pionMinus;
  G4ParticleDefinition* primary2 = proton;

  G4VCollision* theCollision = new G4CollisionMesonBaryon();
  



  G4double mass1 = primary1->GetPDGMass();
  // Formation time
  G4double time1 = 0.;
   // Initial position
  G4double x1 = 0.;
  G4double y1 = 0.;
  G4double z1 = 0.;
 
  G4ThreeVector position1(x1, y1, z1);

  // Select track 2

  G4double mass2 = primary2->GetPDGMass();
  // Formation time
  G4double time2 = 0.;

  // Initial position
  G4double x2 = 0.;
  G4double y2 = 0.;
  G4double z2 = 0.;

  G4ThreeVector position2(x2, y2, z2);


  // Momentum
  //  G4double pzMin = 0.1 * GeV;
  // G4double pzMax = 4.0 * GeV;
  G4double pzMin = 0.527 * GeV;
  G4double pzMax = 4.0 * GeV;


  //    collision->Print();

  G4int n = 195;
  //  G4int n = 195;

  G4double step = (pzMax - pzMin) / n;

  for (G4int i=0; i<n; i++)
    {
      G4double pStep = pzMin + step * i;
      G4ThreeVector p1(0.,0.,pStep);
      G4double e1 = sqrt(mass1*mass1 + p1.mag()*p1.mag());
      G4LorentzVector p41(p1,e1);
      G4KineticTrack trk1(primary1,time1,position1,p41);

      G4ThreeVector p2(0.,0.,0.);
      G4double e2 = sqrt(mass2*mass2 + p2.mag()*p2.mag());
      G4LorentzVector p42(p2,e2);
      G4KineticTrack trk2(primary2,time2,position2,p42);

      G4double sqrts = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();

      // G4double plab = sqrt(sqrts*sqrts - mass1*mass1 - mass2*mass2)/2.0/mass2*mass2;
      // plab = sqrt(plab*plab - mass1*mass1);

      cout << pStep/GeV << "  " << theCollision->CrossSection(trk1, trk2)*GeV << endl;

      
    }


  return EXIT_SUCCESS;
}





