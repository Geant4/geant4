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
#include "G4StringChipsParticleLevelInterface.hh"
#include "globals.hh"
#include "G4Pair.hh"
#include "g4std/list"
#include "g4std/vector"
#include "G4KineticTrackVector.hh"
#include "G4Nucleon.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4LorentzRotation.hh"

G4StringChipsParticleLevelInterface::G4StringChipsParticleLevelInterface()
{
  theEnergyLossPerFermi = 0.7*GeV;
  nop = 164; // ??????
  fractionOfSingleQuasiFreeNucleons = 0.1;
  fractionOfPairedQuasiFreeNucleons = 0.;
  clusteringCoefficient = 2.7;
  temperature = 180.;
  halfTheStrangenessOfSee = 0.3; // = s/d = s/u
  etaToEtaPrime = 0.3;
  fusionToExchange = 100.;
  theInnerCoreDensityCut = 50.;
  
  if(getenv("ChipsParameterTuning"))
  {
    G4cout << "Please enter the energy loss per fermi in GeV"<<G4endl;
    G4cin >> theEnergyLossPerFermi;
    theEnergyLossPerFermi *= GeV;
    G4cout << "Please enter nop"<<G4endl;
    G4cin >> nop;
    G4cout << "Please enter the fractionOfSingleQuasiFreeNucleons"<<G4endl;
    G4cin >> fractionOfSingleQuasiFreeNucleons;
    G4cout << "Please enter the fractionOfPairedQuasiFreeNucleons"<<G4endl;
    G4cin >> fractionOfPairedQuasiFreeNucleons;
    G4cout << "Please enter the clusteringCoefficient"<<G4endl;
    G4cin >> clusteringCoefficient;
    G4cout << "Please enter the temperature"<<G4endl;
    G4cin >> temperature;
    G4cout << "Please enter halfTheStrangenessOfSee"<<G4endl;
    G4cin >> halfTheStrangenessOfSee;
    G4cout << "Please enter the etaToEtaPrime"<<G4endl;
    G4cin >> etaToEtaPrime;
    G4cout << "Please enter the fusionToExchange"<<G4endl;
    G4cin >> fusionToExchange;
    G4cout << "Please enter the cut-off for calculating the nuclear radius in percent"<<G4endl;
    G4cin >> theInnerCoreDensityCut;
  }
}

G4VParticleChange* G4StringChipsParticleLevelInterface::
ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus)
{
  return theModel.ApplyYourself(aTrack, theNucleus);
}

G4ReactionProductVector* G4StringChipsParticleLevelInterface::
Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
{
  // Protection for non physical conditions
  
  if(theSecondaries->size() == 1) 
  {
    G4ReactionProductVector * theFastResult = new G4ReactionProductVector;
    G4ReactionProduct * theFastSec;
    theFastSec = new G4ReactionProduct((*theSecondaries)[0]->GetDefinition());
    G4LorentzVector current4Mom = (*theSecondaries)[0]->Get4Momentum();
    theFastSec->SetTotalEnergy(current4Mom.t());
    theFastSec->SetMomentum(current4Mom.vect());
    theFastResult->push_back(theFastSec);
    return theFastResult;
//    G4Exception("G4StringChipsParticleLevelInterface: Only one particle from String models!");
  }
  
  // target properties needed in constructor of quasmon, and for boosting to
  // target rest frame
  // remove all hit nucleons to get Target code
  theNucleus->StartLoop();
  G4Nucleon * aNucleon;
  G4int resA = 0;
  G4int resZ = 0;
  G4ThreeVector hitMomentum(0,0,0);
  G4double hitMass = 0;
  unsigned int hitCount = 0;
  while((aNucleon = theNucleus->GetNextNucleon()))
  {
    if(!aNucleon->AreYouHit())
    {
      resA++;
      resZ+=G4int (aNucleon->GetDefinition()->GetPDGCharge());
    }
    else
    {
      hitMomentum += aNucleon->GetMomentum().vect();
      hitMass += aNucleon->GetMomentum().m();
      hitCount ++;
    }
  }
  G4int targetPDGCode = 90000000 + 1000*resZ + (resA-resZ);
  G4double targetMass = theNucleus->GetMass();
  targetMass -= hitMass;
  G4double targetEnergy = sqrt(hitMomentum.mag2()+targetMass*targetMass);
  // !! @@ Target should be at rest: hitMomentum=(0,0,0) @@ !! M.K.
  G4LorentzVector targ4Mom(-1.*hitMomentum, targetEnergy);
  
  // Calculate the mean energy lost
  G4Pair<G4double, G4double> theImpact = theNucleus->RefetchImpactXandY();
  G4double impactX = theImpact.first;
  G4double impactY = theImpact.second;
  G4double inpactPar2 = impactX*impactX + impactY*impactY;
  
  G4double radius2 = theNucleus->GetNuclearRadius(theInnerCoreDensityCut*perCent);
  radius2 *= radius2;
  G4double pathlength = 0;
  if(radius2 - inpactPar2>0) pathlength = 2.*sqrt(radius2 - inpactPar2);
  G4double theEnergyLostInFragmentation = theEnergyLossPerFermi*pathlength/fermi;
  
  // now select all particles in range
  G4std::list<G4Pair<G4double, G4KineticTrack *> > theSorted;
  G4std::list<G4Pair<G4double, G4KineticTrack *> >::iterator current;
  for(unsigned int secondary = 0; secondary<theSecondaries->size(); secondary++)
  {
    G4LorentzVector a4Mom = theSecondaries->operator[](secondary)->Get4Momentum();
#ifdef CHIPSdebug
    G4cout <<"ALL STRING particles "<<theSecondaries->operator[](secondary)->GetDefinition()->GetPDGCharge()<<" "
           << theSecondaries->operator[](secondary)->GetDefinition()->GetPDGEncoding()<<" "
	   << a4Mom <<G4endl; 
#endif
    G4double toSort = a4Mom.rapidity();
    G4Pair<G4double, G4KineticTrack *> it;
    it.first = toSort;
    it.second = theSecondaries->operator[](secondary);
    G4bool inserted = false;
    for(current = theSorted.begin(); current!=theSorted.end(); current++)
    {
      if((*current).first > toSort)
      {
	theSorted.insert(current, it);
	inserted = true;
	break;
      }
    }
    if(!inserted)
    {
      theSorted.push_back(it);
    }
  }
  
  G4LorentzVector proj4Mom(0.,0.,0.,0.);
  G4int nD  = 0;
  G4int nU  = 0;
  G4int nS  = 0;
  G4int nAD = 0;
  G4int nAU = 0;
  G4int nAS = 0;
  G4std::list<G4Pair<G4double, G4KineticTrack *> >::iterator firstEscaping = theSorted.begin();
  G4double runningEnergy = 0;
  G4int particleCount = 0;
  G4LorentzVector theLow = (*(theSorted.begin())).second->Get4Momentum();
  G4LorentzVector theHigh;

#ifdef CHIPSdebug
  G4cout << "CHIPS ENERGY LOST "<<theEnergyLostInFragmentation<<G4endl;
  G4cout << "sorted rapidities event start"<<G4endl;
#endif

  G4QHadronVector projHV;
  G4std::vector<G4QContent> theContents;
  G4std::vector<G4LorentzVector *> theMomenta;
  G4ReactionProductVector * theResult = new G4ReactionProductVector;
  G4ReactionProduct * theSec;
  G4KineticTrackVector * secondaries;
  
  for(current = theSorted.begin(); current!=theSorted.end(); current++)
  {
    firstEscaping = current;
    if((*current).second->GetDefinition()->GetQuarkContent(3)!=0 ||
       (*current).second->GetDefinition()->GetAntiQuarkContent(3) !=0)
    {
      G4KineticTrack * aResult = (*current).second;
      G4ParticleDefinition * pdef=aResult->GetDefinition();
      secondaries = NULL;
      if ( pdef->GetPDGWidth() > 0 && pdef->GetPDGLifeTime() < 5E-17*s )
      {
        secondaries = aResult->Decay();
      }
      if ( secondaries == NULL )
      {
        theSec = new G4ReactionProduct(aResult->GetDefinition());
        G4LorentzVector current4Mom = aResult->Get4Momentum();
        current4Mom.boost(targ4Mom.boostVector());
        theSec->SetTotalEnergy(current4Mom.t());
        theSec->SetMomentum(current4Mom.vect());
        theResult->push_back(theSec);
      } 
      else
      {
        for (unsigned int aSecondary=0; aSecondary<secondaries->size(); aSecondary++)
        {
          theSec = new G4ReactionProduct(secondaries->operator[](aSecondary)->GetDefinition());
          G4LorentzVector current4Mom = secondaries->operator[](aSecondary)->Get4Momentum();
          current4Mom.boost(targ4Mom.boostVector());
          theSec->SetTotalEnergy(current4Mom.t());
          theSec->SetMomentum(current4Mom.vect());
          theResult->push_back(theSec);
        }
        G4std::for_each(secondaries->begin(), secondaries->end(), DeleteKineticTrack());
        delete secondaries;
      }
    }

    runningEnergy += (*current).second->Get4Momentum().t();
    if((*current).second->GetDefinition() == G4Proton::Proton())
      runningEnergy-=G4Proton::Proton()->GetPDGMass();
    if((*current).second->GetDefinition() == G4Neutron::Neutron())
      runningEnergy-=G4Neutron::Neutron()->GetPDGMass();

#ifdef CHIPSdebug
      G4cout << "sorted rapidities "<<(*current).second->Get4Momentum().rapidity()<<G4endl;  
#endif

    if(runningEnergy > theEnergyLostInFragmentation) break;
    
#ifdef CHIPSdebug
     G4cout <<"ABSORBED STRING particles "<<(*current).second->GetDefinition()->GetPDGCharge()<<" "
           << (*current).second->GetDefinition()->GetPDGEncoding()<<" "
	   << (*current).second->Get4Momentum() <<G4endl; 
#endif

   // projectile 4-momentum in target rest frame needed in constructor of QHadron
    particleCount++;
    theHigh = (*current).second->Get4Momentum(); 
    proj4Mom = (*current).second->Get4Momentum(); 
    proj4Mom.boost(-1.*targ4Mom.boostVector());  
    nD = (*current).second->GetDefinition()->GetQuarkContent(1);
    nU = (*current).second->GetDefinition()->GetQuarkContent(2);
    nS = (*current).second->GetDefinition()->GetQuarkContent(3);
    nAD = (*current).second->GetDefinition()->GetAntiQuarkContent(1);
    nAU = (*current).second->GetDefinition()->GetAntiQuarkContent(2);
    nAS = (*current).second->GetDefinition()->GetAntiQuarkContent(3);
    G4QContent aProjectile(nD, nU, nS, nAD, nAU, nAS);
//    G4cout << "Quark content: d="<<nD<<", u="<<nU<< ", s="<< nS << G4endl;
//    G4cout << "Anti-quark content: anit-d="<<nAD<<", anti-u="<<nAU<< ", anti-s="<< nAS << G4endl;
//    G4cout << "G4QContent is constructed"<<endl;
    theContents.push_back(aProjectile);
    G4LorentzVector * aVec = new G4LorentzVector(1./MeV*proj4Mom);
    theMomenta.push_back(aVec);
  }
  G4std::vector<G4QContent> theFinalContents;
  G4std::vector<G4LorentzVector*> theFinalMomenta;
  if(theContents.size()<hitCount || 1)
  {
    for(unsigned int hp = 0; hp<theContents.size(); hp++)
    {
      G4QHadron* aHadron = new G4QHadron(theContents[hp], *(theMomenta[hp]) );
      projHV.push_back(aHadron);
    }
  }
  else
  {
    unsigned int hp;
    for(hp=0; hp<hitCount; hp++) 
    {
      G4QContent co(0, 0, 0, 0, 0, 0);
      theFinalContents.push_back(co);
      G4LorentzVector * mo = new G4LorentzVector(0,0,0,0);
      theFinalMomenta.push_back(mo);
    }
    unsigned int running = 0;
    while (running<theContents.size())
    {
      for(hp = 0; hp<hitCount; hp++)
      {
        theFinalContents[hp] +=theContents[running];
	*(theFinalMomenta[hp])+=*(theMomenta[running]);
	running++;
	if(running == theContents.size()) break;
      }
    }
    for(hp = 0; hp<hitCount; hp++)
    {
      G4QHadron* aHadron = new G4QHadron(theFinalContents[hp], *theFinalMomenta[hp]);
      projHV.push_back(aHadron);
    }
  }
  // construct the quasmon
  size_t i;
  for (i=0; i<theFinalMomenta.size(); i++) delete theFinalMomenta[i];
  for (i=0; i<theMomenta.size();      i++) delete theMomenta[i];
  theFinalMomenta.clear();
  theMomenta.clear();

  G4QNucleus::SetParameters(fractionOfSingleQuasiFreeNucleons,
                            fractionOfPairedQuasiFreeNucleons,
			                clusteringCoefficient, 
					fusionToExchange);
  G4Quasmon::SetParameters(temperature,
                           halfTheStrangenessOfSee,
			               etaToEtaPrime);

#ifdef CHIPSdebug
  G4cout << "G4QNucleus parameters "<< fractionOfSingleQuasiFreeNucleons << " "
         << fractionOfPairedQuasiFreeNucleons << " "<< clusteringCoefficient << G4endl;
  G4cout << "G4Quasmon parameters "<< temperature << " "<< halfTheStrangenessOfSee << " "
         << etaToEtaPrime << G4endl;
  G4cout << "The Target PDG code = "<<targetPDGCode<<G4endl;
  G4cout << "The projectile momentum = "<<1./MeV*proj4Mom<<G4endl;
  G4cout << "The target momentum = "<<1./MeV*targ4Mom<<G4endl;
#endif

  // now call chips with this info in place
  G4QHadronVector * output = 0;
  if (particleCount!=0 && resA!=0)
  {
    G4QCHIPSWorld aWorld(nop);              // Create CHIPS World of nop particles
    G4QEnvironment* pan= new G4QEnvironment(projHV, targetPDGCode);
    // clean up particles
    G4std::for_each(projHV.begin(), projHV.end(), DeleteQHadron());
    projHV.clear();
    output = pan->Fragment();
    delete pan;
  }
  else 
  {
    output = new G4QHadronVector;
  }
   
  // Fill the result.
#ifdef CHIPSdebug
  G4cout << "NEXT EVENT"<<endl;
#endif

  // first decay and add all escaping particles.
  for(current = firstEscaping; current!=theSorted.end(); current++)
  {
    G4KineticTrack * aResult = (*current).second;
    G4ParticleDefinition * pdef=aResult->GetDefinition();
    secondaries = NULL;
    if ( pdef->GetPDGWidth() > 0 && pdef->GetPDGLifeTime() < 5E-17*s )
    {
      secondaries = aResult->Decay();
    }
    if ( secondaries == NULL )
    {
      theSec = new G4ReactionProduct(aResult->GetDefinition());
      G4LorentzVector current4Mom = aResult->Get4Momentum();
      current4Mom.boost(targ4Mom.boostVector());
      theSec->SetTotalEnergy(current4Mom.t());
      theSec->SetMomentum(current4Mom.vect());
      theResult->push_back(theSec);
    } 
    else
    {
      for (unsigned int aSecondary=0; aSecondary<secondaries->size(); aSecondary++)
      {
        theSec = new G4ReactionProduct(secondaries->operator[](aSecondary)->GetDefinition());
        G4LorentzVector current4Mom = secondaries->operator[](aSecondary)->Get4Momentum();
        current4Mom.boost(targ4Mom.boostVector());
        theSec->SetTotalEnergy(current4Mom.t());
        theSec->SetMomentum(current4Mom.vect());
        theResult->push_back(theSec);
      }
      G4std::for_each(secondaries->begin(), secondaries->end(), DeleteKineticTrack());
      delete secondaries;
    }
  }
  G4std::for_each(theSecondaries->begin(), theSecondaries->end(), DeleteKineticTrack());
  delete theSecondaries;
    
  // now add the quasmon output
  G4int maxParticle=output->size();
#ifdef CHIPSdebug
  G4cout << "Number of particles from string"<<theResult->size()<<G4endl;
  G4cout << "Number of particles from chips"<<maxParticle<<G4endl;
#endif
  for(G4int particle = 0; particle < maxParticle; particle++)
  {
    if(output->operator[](particle)->GetNFragments() != 0) 
    {
      delete output->operator[](particle);
      continue;
    }
    G4int pdgCode = output->operator[](particle)->GetPDGCode();


#ifdef CHIPSdebug
    G4cerr << "PDG code of chips particle = "<<pdgCode<<G4endl;
#endif

    G4ParticleDefinition * theDefinition;
    // Note that I still have to take care of strange nuclei
    // For this I need the mass calculation, and a changed interface
    // for ion-table ==> work for Hisaya @@@@@@@
    // Then I can sort out the pdgCode. I also need a decau process 
    // for strange nuclei; may be another chips interface
    if(pdgCode>90000000) 
    {
      G4int aZ = (pdgCode-90000000)/1000;
      if (aZ>1000) aZ=aZ%1000;  // patch for strange nuclei, to be repaired @@@@
      G4int anN = pdgCode-90000000-1000*aZ;
      if(anN>1000) anN=anN%1000; // patch for strange nuclei, to be repaired @@@@
      if(pdgCode==91000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(pdgCode==92000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(pdgCode==93000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(pdgCode==94000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(pdgCode==95000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(pdgCode==96000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(pdgCode==97000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(pdgCode==98000000) theDefinition = G4Lambda::LambdaDefinition();
      else if(aZ == 0 && anN == 1) theDefinition = G4Neutron::Neutron();
      else theDefinition = G4ParticleTable::GetParticleTable()->FindIon(aZ,anN+aZ,0,aZ);
    }    
    else theDefinition = G4ParticleTable::GetParticleTable()->FindParticle(pdgCode);

#ifdef CHIPSdebug
    G4cout << "Particle code produced = "<< pdgCode <<G4endl;
#endif

    theSec = new G4ReactionProduct(theDefinition);
    G4LorentzVector current4Mom = output->operator[](particle)->Get4Momentum();
    current4Mom.boost(targ4Mom.boostVector());
    theSec->SetTotalEnergy(current4Mom.t());
    theSec->SetMomentum(current4Mom.vect());
    theResult->push_back(theSec);
    
#ifdef CHIPSdebug
    G4cout <<"CHIPS particles "<<theDefinition->GetPDGCharge()<<" "
           << theDefinition->GetPDGEncoding()<<" "
	   << current4Mom <<G4endl; 
#endif

    delete output->operator[](particle);
  }
  delete output;

#ifdef CHIPSdebug
  G4cout << "Number of particles"<<theResult->size()<<G4endl;
  G4cout << G4endl;
  G4cout << "QUASMON preparation info "
         << 1./MeV*proj4Mom<<" "
	 << 1./MeV*targ4Mom<<" "
	 << nD<<" "<<nU<<" "<<nS<<" "<<nAD<<" "<<nAU<<" "<<nAS<<" "
	 << hitCount<<" "
	 << particleCount<<" "
	 << theLow<<" "
	 << theHigh<<" "
	 << G4endl;
#endif

  return theResult;
} 
