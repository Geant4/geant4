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
#include "G4PartonStringAnnihilator.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4StringInfoDump.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4Proton.hh"
#include "G4ReactionProduct.hh"
#include "G4ReactionProductVector.hh"
#include "G4KineticTrackVector.hh"  
#include "G4AntiProton.hh"

#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4QGSMFragmentation.hh"
#include "G4LundStringFragmentation.hh"
#include "G4PionMinus.hh"

#include "../../string_fragmentation/src/G4VLongitudinalStringDecay.cc"

main()
{
  G4BosonConstructor Bosons;
  Bosons.ConstructParticle();

  G4MesonConstructor Mesons;
  Mesons.ConstructParticle();

  G4LeptonConstructor Leptons;
  Leptons.ConstructParticle();

  G4BaryonConstructor Barions;
  Barions.ConstructParticle();

  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();

  G4PartonStringAnnihilator * theModel = new G4PartonStringAnnihilator;
  G4LundStringFragmentation * aFrag = new G4LundStringFragmentation;
  double fac;
  G4cin >> fac;
//  aFrag->SetMassCut(2.*G4PionMinus::PionMinus()->GetPDGMass()); 
 aFrag->SetSigmaTransverseMomentum(fac*GeV); 
  aFrag->SetVectorMesonProbability(0.9);

  G4ExcitedStringDecay * theDecay = new G4ExcitedStringDecay( aFrag );
//  G4ExcitedStringDecay * theDecay = new G4ExcitedStringDecay( new G4QGSMFragmentation );
  G4StringInfoDump * theDump = new G4StringInfoDump;
  
  G4double pz;
  G4cout << "Please enter the momentum of the anti-proton"<<G4endl;
  G4cin >> pz;
  pz*=GeV;
  int nEvents;
  G4cin >> nEvents;
  G4Fancy3DNucleus * aNucleus = 0;
  G4double aFormationTime = 0;
  G4ThreeVector aPosition(0,0,0);
  G4LorentzVector a4Target(0,0,0,G4Proton::Proton()->GetPDGMass());
  G4LorentzVector a4Projectile(0,0,pz,sqrt(pz*pz + sqr(G4Proton::Proton()->GetPDGMass())));
  G4KineticTrack aTarget(G4Proton::Proton(), aFormationTime, aPosition, a4Target);
  G4KineticTrack aProjectile(G4AntiProton::AntiProton(), aFormationTime, aPosition, a4Projectile);
  
  for(int i=1; i<nEvents; i++)
  {
    G4cerr << "counting events - event #"<<i<<G4endl;
    G4ExcitedString * theString = theModel->GetString(aTarget, aProjectile);
    G4KineticTrackVector * result = theDecay->FragmentString(*theString);
    G4ReactionProductVector* theResult = theDump->Propagate(result, aNucleus);
    
  }
}
