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

#include "G4PreCompoundCascadeInterface.hh"
#include "globals.hh"
#include "G4DynamicParticleVector.hh"
#include "G4IonTable.hh"
#include "G4PreCompoundInuclCollider.hh"
#include "G4IntraNucleiCascader.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4NonEquilibriumEvaporator.hh"
#include "G4BigBanger.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4CollisionOutput.hh"
#include "G4V3DNucleus.hh"
#include "G4Track.hh"
#include "G4Nucleus.hh"
#include "G4NucleiModel.hh"
#include "G4LorentzRotation.hh"


typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;
typedef std::vector<G4InuclNuclei>::iterator nucleiIterator;

G4PreCompoundCascadeInterface::G4PreCompoundCascadeInterface()
  :verboseLevel(0)  {

  if (verboseLevel > 3) {
    G4cout << " >>> G4PreCompoundCascadeInterface::G4PreCompoundCascadeInterface" << G4endl;
  }
}
   
G4ReactionProductVector* G4PreCompoundCascadeInterface::Propagate(G4KineticTrackVector* , 
								  G4V3DNucleus* ) {
  return 0;
}

// #define debug_G4PreCompoundCascadeInterface

G4HadFinalState* G4PreCompoundCascadeInterface::ApplyYourself(const G4HadProjectile& aTrack, 
							      G4Nucleus& theNucleus) {
#ifdef debug_G4PreCompoundCascadeInterface
  static G4int counter(0);
  counter++;
  G4cerr << "Reaction number "<< counter << " "<<aTrack.GetDynamicParticle()->GetDefinition()->GetParticleName()<<" "<< aTrack.GetDynamicParticle()->GetKineticEnergy()<<G4endl;
#endif

  theResult.Clear();

  if (verboseLevel > 3) {
    G4cout << " >>> G4PreCompoundCascadeInterface::ApplyYourself" << G4endl;
  };

  G4double eInit     = 0.0;
  G4double eTot      = 0.0;
  G4double sumBaryon = 0.0;
  G4double sumEnergy = 0.0;

  // Make conversion between native Geant4 and Bertini cascade classes.
  // NOTE: Geant4 units are MeV = 1 and GeV = 1000. Cascade code by default use GeV = 1.

  enum particleType { nuclei      = 0,  proton     = 1,  neutron   = 2,  pionPlus = 3,
                      pionMinus   = 5,  pionZero   = 7,  photon    = 10,
                      kaonPlus    = 11, kaonMinus  = 13, kaonZero  = 15,
                      kaonZeroBar = 17, lambda     = 21, sigmaPlus = 23,
                      sigmaZero   = 25, sigmaMinus = 27, xiZero    = 29, xiMinus  = 31 };

  G4int bulletType = 0;

  // Coding particles 
  if (aTrack.GetDefinition() == G4Proton::Proton()         ) bulletType = proton;
  if (aTrack.GetDefinition() == G4Neutron::Neutron()       ) bulletType = neutron;
  if (aTrack.GetDefinition() == G4PionPlus::PionPlus()     ) bulletType = pionPlus;
  if (aTrack.GetDefinition() == G4PionMinus::PionMinus()   ) bulletType = pionMinus;
  if (aTrack.GetDefinition() == G4PionZero::PionZero()     ) bulletType = pionZero;
  if (aTrack.GetDefinition() == G4Gamma::Gamma()           ) bulletType = photon;
  if (aTrack.GetDefinition() == G4KaonPlus::KaonPlus()     ) bulletType = kaonPlus;
  if (aTrack.GetDefinition() == G4KaonMinus::KaonMinus()   ) bulletType = kaonMinus;
  if (aTrack.GetDefinition() == G4Lambda::Lambda()         ) bulletType = lambda;
  if (aTrack.GetDefinition() == G4SigmaPlus::SigmaPlus()   ) bulletType = sigmaPlus;
  if (aTrack.GetDefinition() == G4SigmaZero::SigmaZero()   ) bulletType = sigmaZero;
  if (aTrack.GetDefinition() == G4SigmaMinus::SigmaMinus() ) bulletType = sigmaMinus;
  if (aTrack.GetDefinition() == G4XiZero::XiZero()         ) bulletType = xiZero;
  if (aTrack.GetDefinition() == G4XiMinus::XiMinus()       ) bulletType = xiMinus;  

  if (aTrack.GetDefinition() == G4KaonZeroLong::KaonZeroLong() ||
      aTrack.GetDefinition() == G4KaonZeroShort::KaonZeroShort() ) {
    if (G4UniformRand() > 0.5) {
      bulletType = kaonZero;
    } else {
      bulletType = kaonZeroBar;
    }
  }

  // Code momentum and energy.
  G4double px,py,pz;
  px=aTrack.Get4Momentum().px() / GeV;
  py=aTrack.Get4Momentum().py() / GeV;
  pz=aTrack.Get4Momentum().pz() / GeV;

  G4LorentzVector projectileMomentum = aTrack.Get4Momentum();
  G4LorentzRotation toZ;
  toZ.rotateZ(-projectileMomentum.phi());
  toZ.rotateY(-projectileMomentum.theta());
  G4LorentzRotation toLabFrame = toZ.inverse();

  std::vector<G4double> momentumBullet(4);
  momentumBullet[0] =0.;
  momentumBullet[1] =0;
  momentumBullet[2] =0;
  momentumBullet[3] =std::sqrt(px*px+py*py+pz*pz);

  G4InuclElementaryParticle *  bullet = new G4InuclElementaryParticle(momentumBullet, bulletType); 

  sumEnergy = bullet->getKineticEnergy(); // In GeV 

  if (bulletType == proton || bulletType == neutron || bulletType == lambda ||
      bulletType == sigmaPlus || bulletType == sigmaZero || bulletType == sigmaMinus ||
      bulletType == xiZero || bulletType == xiMinus) {

    sumBaryon += 1;
  } 

  // Set target
  G4InuclNuclei*   target  = 0;
  G4InuclParticle* targetH = 0;
  // and outcoming particles
  G4DynamicParticle* cascadeParticle = 0;

  std::vector<G4double> targetMomentum(4, 0.0);

  G4double theNucleusA = theNucleus.GetN();

  if ( !(G4int(theNucleusA) == 1) ) {
    target  = new G4InuclNuclei(targetMomentum, 
				theNucleusA, 
				theNucleus.GetZ());
    target->setEnergy();

    std::vector<G4double>  bmom = bullet->getMomentum();
    eInit = std::sqrt(bmom[0] * bmom[0]);
    std::vector<G4double> tmom = target->getMomentum();
    eInit += std::sqrt(tmom[0] * tmom[0]);

    sumBaryon += theNucleusA;

    if (verboseLevel > 2) {
      G4cout << "Bullet:  " << G4endl;  
      bullet->printParticle();
    }
    if (verboseLevel > 2) {
      G4cout << "Target:  " << G4endl;  
      target->printParticle();
    }
  }

  G4CollisionOutput output;

  // Colliders initialisation
  G4ElementaryParticleCollider*   colep = new G4ElementaryParticleCollider;

  G4IntraNucleiCascader*            inc = new G4IntraNucleiCascader; // the actual cascade
  inc->setInteractionCase(1); // Interaction type is particle with nuclei.

  G4NonEquilibriumEvaporator*     noneq = new G4NonEquilibriumEvaporator;
  G4BigBanger*                     bigb = new G4BigBanger;
  G4PreCompoundInuclCollider*  collider = new G4PreCompoundInuclCollider(colep, inc, noneq, bigb);

  G4int  maxTries = 10; // maximum tries for inelastic collision to avoid infinite loop
  G4int  nTries   = 0;  // try counter

  if (G4int(theNucleusA) == 1) { // special treatment for target H(1,1) (proton)

    targetH = new G4InuclElementaryParticle(targetMomentum, 1);

    G4float cutElastic[32];

    cutElastic[proton   ]   = 1.0; // 1 GeV
    cutElastic[neutron  ]   = 1.0;
    cutElastic[lambda]      = 1.0;
    cutElastic[sigmaPlus]   = 1.0;
    cutElastic[sigmaZero]   = 1.0;
    cutElastic[sigmaMinus]  = 1.0;
    cutElastic[xiZero]      = 1.0;
    cutElastic[xiMinus]     = 1.0;

    cutElastic[pionPlus ]   = 0.6; // 0.6 GeV

    cutElastic[kaonPlus ]   = 0.5; // 0.5 GeV
    cutElastic[kaonMinus]   = 0.5;
    cutElastic[kaonMinus]   = 0.5;
    cutElastic[kaonZero]    = 0.5;
    cutElastic[kaonZeroBar] = 0.5;

    cutElastic[pionMinus]   = 0.2; // 0.2 GeV
    cutElastic[pionZero ]   = 0.2;


    if (momentumBullet[3] > cutElastic[bulletType]) { // inelastic collision possible

      do {   // we try to create inelastic interaction
	output = collider->collide(bullet, targetH);
	nTries++;
      } while(
	      (nTries < maxTries) &&
	      (output.getOutgoingParticles().size() == 2 && // elastic: bullet + p = H(1,1) coming out
	       (output.getOutgoingParticles().begin()->type() == bulletType || output.getOutgoingParticles().begin()->type() == proton)
	       )
	      );

    } else { // only elastic collision is energetically possible
      output = collider->collide(bullet, targetH);
    }

    sumBaryon += 1;

    std::vector<G4double> bmom = bullet->getMomentum();
    eInit = std::sqrt(bmom[0] * bmom[0]);
    std::vector<G4double> tmom = targetH->getMomentum();
    eInit += std::sqrt(tmom[0] * tmom[0]);

    if (verboseLevel > 2) {
      G4cout << "Target:  " << G4endl;
      targetH->printParticle();
    }

  } else {  // treat all other targets excepet H(1,1)

    do  // we try to create inelastic interaction
      {
	output = collider->collide(bullet, target );
	nTries++;
      } while (
	       //	      (nTries < maxTries)                                                               &&
	       //(output.getOutgoingParticles().size() + output.getNucleiFragments().size() < 2.5) &&
	       //(output.getOutgoingParticles().size()!=0)                                         &&
	       //(output.getOutgoingParticles().begin()->type()==bullet->type())
	       //);

	       (nTries < maxTries) &&
	       output.getOutgoingParticles().size() == 1 &&     // we retry when elastic collision happened
               output.getNucleiFragments().size() == 1 &&            
	       output.getOutgoingParticles().begin()->type() == bullet->type() &&
	       output.getNucleiFragments().begin()->getA() == target->getA() && 
	       output.getNucleiFragments().begin()->getZ() == target->getZ() 
	       );
  }
 
  if (verboseLevel > 1) 
    {
      G4cout << " Cascade output: " << G4endl;
      output.printCollisionOutput();
    }
  
  // Convert cascade data to use hadronics interface
  std::vector<G4InuclNuclei>             nucleiFragments = output.getNucleiFragments();
  std::vector<G4InuclElementaryParticle> particles =       output.getOutgoingParticles();

  theResult.SetStatusChange(stopAndKill);

  if (!particles.empty()) { 
    particleIterator ipart;
    G4int outgoingParticle;

    for (ipart = particles.begin(); ipart != particles.end(); ipart++) {
      outgoingParticle = ipart->type();
      std::vector<G4double> mom = ipart->getMomentum();
      eTot   += std::sqrt(mom[0] * mom[0]);

      G4double ekin = ipart->getKineticEnergy() * GeV;
      G4ThreeVector aMom(mom[1], mom[2], mom[3]);
      aMom = aMom.unit();

      if (ipart->baryon() ) {
	sumBaryon -= 1;
      }

      sumEnergy -= ekin / GeV;

      switch(outgoingParticle) {

      case proton: 
#ifdef debug_G4PreCompoundCascadeInterface
	G4cerr << "proton " << counter << " " << aMom << " " << ekin << G4endl;
#endif
	cascadeParticle = new G4DynamicParticle(G4Proton::ProtonDefinition(), aMom, ekin);
	break; 

      case neutron: 

#ifdef debug_G4PreCompoundCascadeInterface
	G4cerr << "neutron "<< counter<<" "<<aMom<<" "<<  ekin<<G4endl;
#endif
	cascadeParticle = new G4DynamicParticle(G4Neutron::NeutronDefinition(), aMom, ekin);
	break;

      case pionPlus: 
	cascadeParticle = new G4DynamicParticle(G4PionPlus::PionPlusDefinition(), aMom, ekin);

#ifdef debug_G4PreCompoundCascadeInterface
	G4cerr << "pionPlus "<< counter<<" "<<aMom<<" "<<  ekin<<G4endl;
#endif
	break;

      case pionMinus:
	cascadeParticle = new G4DynamicParticle(G4PionMinus::PionMinusDefinition(), aMom, ekin);

#ifdef debug_G4PreCompoundCascadeInterface
	G4cerr << "pionMinus "<< counter<<" "<<aMom<<" "<<  ekin<<G4endl;
#endif
	break;

      case pionZero: 
	cascadeParticle = new G4DynamicParticle(G4PionZero::PionZeroDefinition(), aMom, ekin);

#ifdef debug_G4PreCompoundCascadeInterface
	G4cerr << "pionZero "<< counter<<" "<<aMom<<" "<<  ekin<<G4endl;
#endif
	break;

      case photon: 
	cascadeParticle = new G4DynamicParticle(G4Gamma::Gamma(), aMom, ekin);

#ifdef debug_G4PreCompoundCascadeInterface
	G4cerr << "photon "<< counter<<" "<<aMom<<" "<<  ekin<<G4endl;
#endif
	break;


      case kaonPlus:
        cascadeParticle = new G4DynamicParticle(G4KaonPlus::KaonPlusDefinition(), aMom, ekin);
        break;

      case kaonMinus:
        cascadeParticle = new G4DynamicParticle(G4KaonMinus::KaonMinusDefinition(), aMom, ekin);
        break;

      case kaonZero:
        if (G4UniformRand() > 0.5) {
          cascadeParticle = new G4DynamicParticle(G4KaonZeroLong::KaonZeroLongDefinition(), aMom, ekin);
        } else {
          cascadeParticle = new G4DynamicParticle(G4KaonZeroShort::KaonZeroShortDefinition(), aMom, ekin);
        }
        break;

      case kaonZeroBar:
        if (G4UniformRand() > 0.5) {
          cascadeParticle = new G4DynamicParticle(G4KaonZeroLong::KaonZeroLongDefinition(), aMom, ekin);
        } else {
          cascadeParticle = new G4DynamicParticle(G4KaonZeroShort::KaonZeroShortDefinition(), aMom, ekin);
        }
        break;

      case lambda:
        cascadeParticle = new G4DynamicParticle(G4Lambda::LambdaDefinition(), aMom, ekin);
        break;

      case sigmaPlus:
        cascadeParticle = new G4DynamicParticle(G4SigmaPlus::SigmaPlusDefinition(), aMom, ekin);
        break;

      case sigmaZero:
        cascadeParticle = new G4DynamicParticle(G4SigmaZero::SigmaZeroDefinition(), aMom, ekin);
        break;

      case sigmaMinus:
        cascadeParticle = new G4DynamicParticle(G4SigmaMinus::SigmaMinusDefinition(), aMom, ekin);
        break;

      case xiZero:
        cascadeParticle = new G4DynamicParticle(G4XiZero::XiZeroDefinition(), aMom, ekin);
        break;

      case xiMinus:
        cascadeParticle = new G4DynamicParticle(G4XiMinus::XiMinusDefinition(), aMom, ekin);
        break;

      default:
        G4cout << " ERROR: G4PreCompoundCascadeInterface::Propagate undefined particle type" << G4endl;
      }

      cascadeParticle->Set4Momentum(cascadeParticle->Get4Momentum()*=toLabFrame);
      theResult.AddSecondary(cascadeParticle); 
    }
  }

  // get nuclei fragments
  G4DynamicParticle * aFragment = 0;
  G4ParticleDefinition * aIonDef = 0;
  G4ParticleTable *theTableOfParticles = G4ParticleTable::GetParticleTable();

  if (!nucleiFragments.empty()) { 
    nucleiIterator ifrag;

    for (ifrag = nucleiFragments.begin(); ifrag != nucleiFragments.end(); ifrag++) 
      {
	G4double eKin = ifrag->getKineticEnergy() * GeV;
	std::vector<G4double> mom = ifrag->getMomentum();
        eTot   += std::sqrt(mom[0] * mom[0]);

	G4ThreeVector aMom(mom[1], mom[2], mom[3]);
	aMom = aMom.unit();

	// hpw @@@ ==> Should be zero: G4double fragmentExitation = ifrag->getExitationEnergyInGeV();

	if (verboseLevel > 2) {
	  G4cout << " Nuclei fragment: " << G4endl;
	  ifrag->printParticle();
	}

	G4int A = G4int(ifrag->getA());
	G4int Z = G4int(ifrag->getZ());
	aIonDef = theTableOfParticles->FindIon(Z, A, 0, Z);
      
	aFragment =  new G4DynamicParticle(aIonDef, aMom, eKin);

	sumBaryon -= A;
	sumEnergy -= eKin / GeV;

        aFragment->Set4Momentum(aFragment->Get4Momentum()*=toLabFrame);
	theResult.AddSecondary(aFragment); 
      }
  }

  if (verboseLevel > 2) {
    if (sumBaryon != 0) {
      G4cout << "ERROR: no baryon number conservation, sum of baryons = " << sumBaryon << G4endl;
    }

    if (sumEnergy > 0.01 ) {
      G4cout << "Kinetic energy conservation violated by " << sumEnergy << " GeV" << G4endl;
    }
     
    G4cout << "Total energy conservation at level ~" << (eInit - eTot) * GeV << " MeV" << G4endl;
    
    if (sumEnergy < -5.0e-5 ) { // 0.05 MeV
      G4cout << "FATAL ERROR: energy created  " << sumEnergy * GeV << " MeV" << G4endl;
    }
  }

  delete bullet;
  delete colep;
  delete inc;
  delete noneq; 
  delete bigb;
  delete collider;

  if(target != 0) delete target;
  if(targetH != 0) delete targetH;
  // if(cascadeParticle != 0) delete cascadeParticle;
  // if(aFragment != 0) delete aFragment;

  return &theResult;
  }
