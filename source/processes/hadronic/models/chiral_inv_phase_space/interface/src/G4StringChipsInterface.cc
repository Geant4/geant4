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
// !!! Was used in QBBC PL, NOW it is not. Must be absolete !!!
// ============================================================

//#define debug
//#define pdebug

#include "G4StringChipsInterface.hh"
#include "globals.hh"
#include <utility>
#include <list>
#include "G4KineticTrackVector.hh"
#include "G4Nucleon.hh"
#include "G4LorentzRotation.hh"

G4StringChipsInterface::G4StringChipsInterface()
{
#ifdef CHIPSdebug
  G4cout << "Please enter the energy loss per fermi in GeV"<<G4endl;
  G4cin >> theEnergyLossPerFermi;
#endif
#ifdef debug
  G4cout<<"G4StringChipsInterface::Constructor is called"<<G4endl;
#endif
  theEnergyLossPerFermi = 0.5*GeV;
  // theEnergyLossPerFermi = 1.*GeV;
}

G4HadFinalState* G4StringChipsInterface::
ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& theNucleus)
{
#ifdef debug
  G4cout<<"G4StringChipsInterface::ApplyYourself is called"<<G4endl;
#endif
  return theModel.ApplyYourself(aTrack, theNucleus);
}

G4ReactionProductVector* G4StringChipsInterface::
Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
{
#ifdef debug
  G4cout<<"G4StringChipsInterface::Propagate is called"<<G4endl;
#endif
  // Protection for non physical conditions
  
  if(theSecondaries->size() == 1) 
    throw G4HadronicException(__FILE__, __LINE__, "G4StringChipsInterface: Only one particle from String models!");
  
  // Calculate the mean energy lost
  std::pair<G4double, G4double> theImpact = theNucleus->RefetchImpactXandY();
  G4double impactX = theImpact.first;
  G4double impactY = theImpact.second;
  G4double inpactPar2 = impactX*impactX + impactY*impactY;
  
  G4double radius2 = theNucleus->GetNuclearRadius(5*perCent);
  radius2 *= radius2;
  G4double pathlength = 0;
  if(radius2 - inpactPar2>0) pathlength = 2.*std::sqrt(radius2 - inpactPar2);
  G4double theEnergyLostInFragmentation = theEnergyLossPerFermi*pathlength/fermi;
  
  // now select all particles in range
  std::list<std::pair<G4double, G4KineticTrack *> > theSorted;
  std::list<std::pair<G4double, G4KineticTrack *> >::iterator current;
  for(unsigned int secondary = 0; secondary<theSecondaries->size(); secondary++)
  {
    G4LorentzVector a4Mom = theSecondaries->operator[](secondary)->Get4Momentum();

#ifdef CHIPSdebug
    G4cout <<"ALL STRING particles "<<theSecondaries->operator[](secondary)->GetDefinition()->GetPDGCharge()<<" "
           << theSecondaries->operator[](secondary)->GetDefinition()->GetPDGEncoding()<<" "
	   << a4Mom <<G4endl; 
#endif

    G4double toSort = a4Mom.rapidity();
    std::pair<G4double, G4KineticTrack *> it;
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
      theSorted.push_front(it);
    }
  }
  
  G4LorentzVector proj4Mom(0.,0.,0.,0.);
  // @@ Use the G4QContent class, which is exactly (nD,nU,nS,nAD,nAU,nAS) !
  // The G4QContent class is a basic clas in CHIPS (not PDG Code as in GEANT4),
  // so in CHIPS on can has a hadronic obgect (Quasmon) with any Quark Content.
  // As a simple extantion for the hadron (which is a special case for Hadron)
  // there is a clas G4QChipolino, which is a Quasmon, which can decay in two
  // hadrons. In future the three-hadron class can be added...
  G4int nD  = 0;
  G4int nU  = 0;
  G4int nS  = 0;
  G4int nAD = 0;
  G4int nAU = 0;
  G4int nAS = 0;
  std::list<std::pair<G4double, G4KineticTrack *> >::iterator firstEscaping = theSorted.begin();
  G4double runningEnergy = 0;
  G4int particleCount = 0;
  G4LorentzVector theLow = (*(theSorted.begin())).second->Get4Momentum();
  G4LorentzVector theHigh;

#ifdef CHIPSdebug
  G4cout << "CHIPS ENERGY LOST "<<theEnergyLostInFragmentation<<G4endl;
  G4cout << "sorted rapidities event start"<<G4endl;
#endif 

  G4ReactionProductVector * theResult = new G4ReactionProductVector;
  G4ReactionProduct * theSec;
  G4KineticTrackVector * secondaries;
#ifdef pdebug
  G4cout<<"G4StringChipsInterface::Propagate: Absorption"<<G4endl; 
#endif
  
  // first decay and add all escaping particles.
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
          theSec->SetTotalEnergy(current4Mom.t());
          theSec->SetMomentum(current4Mom.vect());
          theResult->push_back(theSec);
        }
        std::for_each(secondaries->begin(), secondaries->end(), DeleteKineticTrack());
        delete secondaries;
      }
      continue;
    }
    runningEnergy += (*current).second->Get4Momentum().t();
    if(runningEnergy > theEnergyLostInFragmentation) break;
    
#ifdef CHIPSdebug
    G4cout <<"ABSORBED STRING particles "<<current->second->GetDefinition()->GetPDGCharge()<<" "
           << current->second->GetDefinition()->GetPDGEncoding()<<" "
	          << current->second->Get4Momentum() <<G4endl; 
#endif
#ifdef pdebug
    G4cout<<"G4StringChipsInterface::Propagate:C="
          <<current->second->GetDefinition()->GetPDGCharge()<<", PDG="
          << current->second->GetDefinition()->GetPDGEncoding()<<", 4M="
	         << current->second->Get4Momentum() <<G4endl; 
#endif

    // projectile 4-momentum needed in constructor of quasmon
    particleCount++;
    theHigh = (*current).second->Get4Momentum(); 
    proj4Mom += (*current).second->Get4Momentum(); 

#ifdef CHIPSdebug
    G4cout << "sorted rapidities "<<current->second->Get4Momentum().rapidity()<<G4endl;
#endif
    
     // projectile quark contents needed for G4QContent construction (@@ ? -> G4QContent class)
    nD  += (*current).second->GetDefinition()->GetQuarkContent(1);
    nU  += (*current).second->GetDefinition()->GetQuarkContent(2);
    nS  += (*current).second->GetDefinition()->GetQuarkContent(3);
    nAD += (*current).second->GetDefinition()->GetAntiQuarkContent(1);
    nAU += (*current).second->GetDefinition()->GetAntiQuarkContent(2);
    nAS += (*current).second->GetDefinition()->GetAntiQuarkContent(3);
  }
  // construct G4QContent

#ifdef CHIPSdebug
  G4cout << "Quark content: d="<<nD<<", u="<<nU<< ", s="<< nS << G4endl;
  G4cout << "Anti-quark content: anit-d="<<nAD<<", anti-u="<<nAU<< ", anti-s="<< nAS << G4endl;
#endif

  G4QContent theProjectiles(nD, nU, nS, nAD, nAU, nAS);

#ifdef CHIPSdebug
  G4cout << "G4QContent is constructed"<<endl;
#endif
  
  // target properties needed in constructor of quasmon
  // remove all hit nucleons to get Target code
  theNucleus->StartLoop();
  G4Nucleon * aNucleon;
  G4int resA = 0;
  G4int resZ = 0;
  G4ThreeVector hitMomentum(0,0,0);
  G4double hitMass = 0;
  G4int hitCount = 0;
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
  G4double targetEnergy = std::sqrt(hitMomentum.mag2()+targetMass*targetMass);
  // !! @@ Target should be at rest: hitMomentum=(0,0,0) @@ !! M.K.
  G4LorentzVector targ4Mom(-1.*hitMomentum, targetEnergy);
  
  // construct the quasmon
  //  G4int nop = 122; // clusters up to Alpha cluster
  G4int nop = 152; // V.Ivanchenko set the same parameter as for all CHIPS models
  G4double fractionOfSingleQuasiFreeNucleons = 0.5; // It is A-dependent (C=.85, U=.40)
  G4double fractionOfPairedQuasiFreeNucleons = 0.05;
  G4double clusteringCoefficient = 5.;
  G4double temperature = 180.;
  G4double halfTheStrangenessOfSee = 0.3; // = s/d = s/u
  G4double etaToEtaPrime = 0.3;

  G4QNucleus::SetParameters(fractionOfSingleQuasiFreeNucleons,
                            fractionOfPairedQuasiFreeNucleons,
			    clusteringCoefficient);
  G4Quasmon::SetParameters(temperature,
                           halfTheStrangenessOfSee,
			   etaToEtaPrime);

#ifdef CHIPSdebug
  G4cout << "G4QNucleus parameters "<< fractionOfSingleQuasiFreeNucleons << " "
         << fractionOfPairedQuasiFreeNucleons << " "<< clusteringCoefficient << G4endl;
  G4cout << "G4Quasmon parameters "<< temperature << " "<< halfTheStrangenessOfSee << " "
         <<etaToEtaPrime << G4endl;
  G4cout << "The Target PDG code = "<<targetPDGCode<<G4endl;
  G4cout << "The projectile momentum = "<<1./MeV*proj4Mom<<G4endl;
  G4cout << "The target momentum = "<<1./MeV*targ4Mom<<G4endl;
#endif

  
  // Chips expects all in target rest frame, along z.
  // G4QCHIPSWorld aWorld(nop);              // Create CHIPS World of nop particles
  G4QCHIPSWorld::Get()->GetParticles(nop);
  G4QHadronVector projHV;
  // target rest frame
  proj4Mom.boost(-1.*targ4Mom.boostVector());
  // now go along z
  G4LorentzRotation toZ;
  toZ.rotateZ(-1*proj4Mom.phi());
  toZ.rotateY(-1*proj4Mom.theta());
  proj4Mom = toZ*proj4Mom;
  G4LorentzRotation toLab(toZ.inverse());

#ifdef CHIPSdebug
  G4cout << "a Boosted projectile vector along z"<<proj4Mom<<" "<<proj4Mom.mag()<<G4endl;
#endif
  
  G4QHadron* iH = new G4QHadron(theProjectiles, 1./MeV*proj4Mom);
  projHV.push_back(iH);

  // now call chips with this info in place
  G4QHadronVector * output = 0;
  if (particleCount!=0)
  {
    G4QEnvironment* pan= new G4QEnvironment(projHV, targetPDGCode);
    try
    {
      output = pan->Fragment();
    }
    catch(G4HadronicException & aR)
    {
      G4cerr << "Exception thrown passing through G4ChiralInvariantPhaseSpace "<<G4endl;
      G4cerr << " targetPDGCode = "<< targetPDGCode <<G4endl;
      G4cerr << " Dumping the information in the pojectile list"<<G4endl;
      for(size_t i=0; i< projHV.size(); i++)
      {
	       G4cerr <<"  Incoming 4-momentum and PDG code of "<<i<<"'th hadron: "
               <<" "<< projHV[i]->Get4Momentum()<<" "<<projHV[i]->GetPDGCode()<<G4endl;
      }
      throw;
    }
    std::for_each(projHV.begin(), projHV.end(), DeleteQHadron());
    projHV.clear();
    delete pan;
  }
  else output = new G4QHadronVector;
   
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
      current4Mom = toLab*current4Mom;
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
        current4Mom = toLab*current4Mom;
        current4Mom.boost(targ4Mom.boostVector());
        theSec->SetTotalEnergy(current4Mom.t());
        theSec->SetMomentum(current4Mom.vect());
        theResult->push_back(theSec);
      }
      std::for_each(secondaries->begin(), secondaries->end(), DeleteKineticTrack());
      delete secondaries;
    }
  }
  std::for_each(theSecondaries->begin(), theSecondaries->end(), DeleteKineticTrack());
  delete theSecondaries;
    
  // now add the quasmon output
#ifdef CHIPSdebug
  G4cout << "Number of particles from string"<<theResult->size()<<G4endl;
  G4cout << "Number of particles from chips"<<output->size()<<G4endl;
#endif

  for(unsigned int particle = 0; particle < output->size(); particle++)
  {
    if(output->operator[](particle)->GetNFragments() != 0) 
    {
      delete output->operator[](particle);
      continue;
    }
    theSec = new G4ReactionProduct;  
    G4int pdgCode = output->operator[](particle)->GetPDGCode();
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
      theDefinition = G4ParticleTable::GetParticleTable()->FindParticle(output->operator[](particle)->GetPDGCode());
    }
#ifdef CHIPSdebug
    G4cout << "Particle code produced = "<< pdgCode <<G4endl;
#endif

    theSec = new G4ReactionProduct(theDefinition);
    G4LorentzVector current4Mom = output->operator[](particle)->Get4Momentum();
    current4Mom = toLab*current4Mom;
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
  // @@ G4QContent has even the out option!
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
