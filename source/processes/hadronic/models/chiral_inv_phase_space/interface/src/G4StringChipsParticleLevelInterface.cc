#include "G4StringChipsParticleLevelInterface.hh"
#include "globals.hh"
#include "G4Pair.hh"
#include "g4std/list"
#include "G4KineticTrackVector.hh"
#include "G4Nucleon.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4LorentzRotation.hh"

G4StringChipsParticleLevelInterface::G4StringChipsParticleLevelInterface()
{
  G4cout << "Please enter the energy loss per fermi in GeV"<<G4endl;
  G4cin >> theEnergyLossPerFermi;
  theEnergyLossPerFermi *= GeV;
  // theEnergyLossPerFermi = 1.*GeV;
  G4cout << "Please enter nop"<<G4endl;
  G4cin >> nop;
  // G4int nop = 223; // ??????
  G4cout << "Please enter the fractionOfSingleQuasiFreeNucleons"<<G4endl;
  G4cin >> fractionOfSingleQuasiFreeNucleons;
  // G4double fractionOfSingleQuasiFreeNucleons = 0.45;
  G4cout << "Please enter the fractionOfPairedQuasiFreeNucleons"<<G4endl;
  G4cin >> fractionOfPairedQuasiFreeNucleons;
  // G4double fractionOfPairedQuasiFreeNucleons = 0.15;
  G4cout << "Please enter the clusteringCoefficient"<<G4endl;
  G4cin >> clusteringCoefficient;
  //G4double clusteringCoefficient = 5.;
  G4cout << "Please enter the temperature"<<G4endl;
  G4cin >> temperature;
  G4double temperature = 180.;
  G4cout << "Please enter halfTheStrangenessOfSee"<<G4endl;
  G4cin >> halfTheStrangenessOfSee;
  //G4double halfTheStrangenessOfSee = 0.1; // = s/d = s/u
  G4cout << "Please enter the etaToEtaPrime"<<G4endl;
  G4cin >> etaToEtaPrime;
  //G4double etaToEtaPrime = 0.3;
  G4cout << "Please enter the fusionToExchange"<<G4endl;
  G4cin >> fusionToExchange;
  //G4double fusionToExchange = 0.000001;
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
  
  if(theSecondaries->length() == 1) 
    G4Exception("G4StringChipsParticleLevelInterface: Only one particle from String models!");
  
  // target properties needed in constructor of quasmon, and for boosting to
  // target rest frame
  // remove all hit nucleons to get Target code
  theNucleus->StartLoop();
  G4Nucleon * aNucleon;
  G4int resA = 0;
  G4int resZ = 0;
  G4ThreeVector hitMomentum(0,0,0);
  G4double hitMass = 0;
  G4int hitCount = 0;
  while(aNucleon = theNucleus->GetNextNucleon())
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
  
  G4double radius2 = theNucleus->GetNuclearRadius(5*perCent);
  radius2 *= radius2;
  G4double pathlength = 0;
  if(radius2 - inpactPar2>0) pathlength = 2.*sqrt(radius2 - inpactPar2);
  G4double theEnergyLostInFragmentation = theEnergyLossPerFermi*pathlength/fermi;
  
  // now select all particles in range
  G4std::list<G4Pair<G4double, G4KineticTrack *> > theSorted;
  G4std::list<G4Pair<G4double, G4KineticTrack *> >::iterator current;
  for(G4int secondary = 0; secondary<theSecondaries->length(); secondary++)
  {
    G4LorentzVector a4Mom = theSecondaries->at(secondary)->Get4Momentum();
    G4cout <<"ALL STRING particles "<<theSecondaries->at(secondary)->GetDefinition()->GetPDGCharge()<<" "
           << theSecondaries->at(secondary)->GetDefinition()->GetPDGEncoding()<<" "
	   << a4Mom <<G4endl; 
    G4double toSort = a4Mom.rapidity();
    G4Pair<G4double, G4KineticTrack *> it;
    it.first = toSort;
    it.second = theSecondaries->at(secondary);
    G4bool inserted = false;
    for(current = theSorted.begin(); current!=theSorted.end(); current++)
    {
      if(current->first > toSort)
      {
	theSorted.insert(current, it);
	inserted = true;
	break;
      }
    }
    if(!inserted)
    {
      theSorted.push_front(it);
    }
  }
  
  G4LorentzVector proj4Mom = (0,0,0,0);
  G4int nD  = 0;
  G4int nU  = 0;
  G4int nS  = 0;
  G4int nAD = 0;
  G4int nAU = 0;
  G4int nAS = 0;
  G4std::list<G4Pair<G4double, G4KineticTrack *> >::iterator firstEscaping = theSorted.begin();
  G4double runningEnergy = 0;
  G4int particleCount = 0;
  G4LorentzVector theLow = theSorted.begin()->second->Get4Momentum();
  G4LorentzVector theHigh;
  G4cout << "CHIPS ENERGY LOST "<<theEnergyLostInFragmentation<<G4endl;
  G4cout << "sorted rapidities event start"<<G4endl;
  G4QHadronVector projHV;
  G4std::vector<G4QContent> theContents;
  G4std::vector<G4LorentzVector> theMomenta;
  for(current = theSorted.begin(); current!=theSorted.end(); current++)
  {
    firstEscaping = current;
    runningEnergy += current->second->Get4Momentum().t();
    if(current->second->GetDefinition() == G4Proton::Proton())
      runningEnergy-=G4Proton::Proton()->GetPDGMass();
    if(current->second->GetDefinition() == G4Neutron::Neutron())
      runningEnergy-=G4Neutron::Neutron()->GetPDGMass();
    if(runningEnergy > theEnergyLostInFragmentation) break;
    
     G4cout <<"ABSORBED STRING particles "<<current->second->GetDefinition()->GetPDGCharge()<<" "
           << current->second->GetDefinition()->GetPDGEncoding()<<" "
	   << current->second->Get4Momentum() <<G4endl; 
   // projectile 4-momentum in target rest frame needed in constructor of QHadron
    particleCount++;
    theHigh = current->second->Get4Momentum(); 
    proj4Mom = current->second->Get4Momentum(); 
    proj4Mom.boost(-1.*targ4Mom.boostVector());  
//     G4cout << "sorted rapidities "<<current->second->Get4Momentum().rapidity()<<G4endl;  
    nD = current->second->GetDefinition()->GetQuarkContent(1);
    nU = current->second->GetDefinition()->GetQuarkContent(2);
    nS = current->second->GetDefinition()->GetQuarkContent(3);
    nAD = current->second->GetDefinition()->GetAntiQuarkContent(1);
    nAU = current->second->GetDefinition()->GetAntiQuarkContent(2);
    nAS = current->second->GetDefinition()->GetAntiQuarkContent(3);
    G4QContent aProjectile(nD, nU, nS, nAD, nAU, nAS);
//    G4cout << "Quark content: d="<<nD<<", u="<<nU<< ", s="<< nS << G4endl;
//    G4cout << "Anti-quark content: anit-d="<<nAD<<", anti-u="<<nAU<< ", anti-s="<< nAS << G4endl;
//    G4cout << "G4QContent is constructed"<<endl;
    theContents.push_back(aProjectile);
    theMomenta.push_back(1./MeV*proj4Mom);
  }
  G4std::vector<G4QContent> theFinalContents;
  G4std::vector<G4LorentzVector> theFinalMomenta;
  if(theContents.size()<hitCount || 1)
  {
    for(G4int hp = 0; hp<theContents.size(); hp++)
    {
      G4QHadron* aHadron = new G4QHadron(theContents[hp], theMomenta[hp]);
      projHV.insert(aHadron);
    }
  }
  else
  {
    G4int hp;
    for(hp=0; hp<hitCount; hp++) 
    {
      G4QContent co(0, 0, 0, 0, 0, 0);
      theFinalContents.push_back(co);
      G4LorentzVector mo(0,0,0,0);
      theFinalMomenta.push_back(mo);
    }
    G4int running = 0;
    while (running<theContents.size())
    {
      for(hp = 0; hp<hitCount; hp++)
      {
        theFinalContents[hp] +=theContents[running];
	theFinalMomenta[hp]+=theMomenta[running];
	running++;
	if(running == theContents.size()) break;
      }
    }
    for(hp = 0; hp<hitCount; hp++)
    {
      G4QHadron* aHadron = new G4QHadron(theFinalContents[hp], theFinalMomenta[hp]);
      projHV.insert(aHadron);
    }
  }
  // construct the quasmon

  G4QNucleus::SetParameters(fractionOfSingleQuasiFreeNucleons,
                            fractionOfPairedQuasiFreeNucleons,
			                clusteringCoefficient, 
					fusionToExchange);
  G4Quasmon::SetParameters(temperature,
                           halfTheStrangenessOfSee,
			               etaToEtaPrime);
  G4cout << "G4QNucleus parameters "<< fractionOfSingleQuasiFreeNucleons << " "
         << fractionOfPairedQuasiFreeNucleons << " "<< clusteringCoefficient << G4endl;
  G4cout << "G4Quasmon parameters "<< temperature << " "<< halfTheStrangenessOfSee << " "
         << etaToEtaPrime << G4endl;
  G4cout << "The Target PDG code = "<<targetPDGCode<<G4endl;
  G4cout << "The projectile momentum = "<<1./MeV*proj4Mom<<G4endl;
  G4cout << "The target momentum = "<<1./MeV*targ4Mom<<G4endl;

  // now call chips with this info in place
  G4QHadronVector * output = 0;
  if (particleCount!=0)
  {
    G4QCHIPSWorld aWorld(nop);              // Create CHIPS World of nop particles
    G4QEnvironment* pan= new G4QEnvironment(projHV, targetPDGCode);
    // clean up particles
    projHV.clearAndDestroy();
    output = pan->Fragment();
    delete pan;
  }
  else 
  {
    output = new G4QHadronVector;
  }
   
  // Fill the result.
  G4ReactionProductVector * theResult = new G4ReactionProductVector;
  G4ReactionProduct * theSec;
  G4cout << "NEXT EVENT"<<endl;
  
  // first decay and add all escaping particles.
  G4KineticTrackVector * secondaries;
  for(current = firstEscaping; current!=theSorted.end(); current++)
  {
    G4KineticTrack * aResult = current->second;
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
      theResult->insert(theSec);
    } 
    else
    {
      for (G4int aSecondary=0; aSecondary<secondaries->entries(); aSecondary++)
      {
        theSec = new G4ReactionProduct(secondaries->at(aSecondary)->GetDefinition());
        G4LorentzVector current4Mom = secondaries->at(aSecondary)->Get4Momentum();
        current4Mom.boost(targ4Mom.boostVector());
        theSec->SetTotalEnergy(current4Mom.t());
        theSec->SetMomentum(current4Mom.vect());
        theResult->insert(theSec);
      }
      secondaries->clearAndDestroy();
      delete secondaries;
    }
  }
  theSecondaries->clearAndDestroy();
  delete theSecondaries;
    
  // now add the quasmon output
  G4cout << "Number of particles from string"<<theResult->length()<<G4endl;
  G4cout << "Number of particles from chips"<<output->length()<<G4endl;
  for(G4int particle = 0; particle < output->length(); particle++)
  {
    if(output->at(particle)->GetNFragments() != 0) 
    {
      delete output->at(particle);
      continue;
    }
    theSec = new G4ReactionProduct;  
    G4int pdgCode = output->at(particle)->GetPDGCode();
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
    else
    {
      theDefinition = G4ParticleTable::GetParticleTable()->FindParticle(output->at(particle)->GetPDGCode());
    }
    G4cout << "Particle code produced = "<< pdgCode <<G4endl;
    theSec = new G4ReactionProduct(theDefinition);
    G4LorentzVector current4Mom = output->at(particle)->Get4Momentum();
    current4Mom.boost(targ4Mom.boostVector());
    theSec->SetTotalEnergy(current4Mom.t());
    theSec->SetMomentum(current4Mom.vect());
    theResult->insert(theSec);
    
    G4cout <<"CHIPS particles "<<theDefinition->GetPDGCharge()<<" "
           << theDefinition->GetPDGEncoding()<<" "
	   << current4Mom <<G4endl; 
    delete output->at(particle);
  }
  delete output;
  G4cout << "Number of particles"<<theResult->length()<<G4endl;
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
  return theResult;
} 
