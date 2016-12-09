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
#include <utility>

#include "G4QGSParticipants.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4LorentzVector.hh"

// Class G4QGSParticipants 

// HPW Feb 1999
// Promoting model parameters from local variables class properties
G4ThreadLocal G4int G4QGSParticipants_NPart = 0;

G4QGSParticipants::G4QGSParticipants() : theDiffExcitaton(), //0.7*GeV, 250*MeV, 250*MeV),
		ModelMode(SOFT),
		//nCutMax(7),ThresholdParameter(0.45*GeV),
		nCutMax(7),ThresholdParameter(0.000*GeV),
		QGSMThreshold(3*GeV),theNucleonRadius(1.5*fermi)
{
}

G4QGSParticipants::G4QGSParticipants(const G4QGSParticipants &right)
: G4VParticipants(),ModelMode(right.ModelMode), nCutMax(right.nCutMax),
  ThresholdParameter(right.ThresholdParameter), QGSMThreshold(right.QGSMThreshold),
  theNucleonRadius(right.theNucleonRadius)
{
}


G4QGSParticipants::~G4QGSParticipants()
{
}

void G4QGSParticipants::BuildInteractions(const G4ReactionProduct  &thePrimary) 
{
  // Find the collisions and collition conditions
  G4VSplitableHadron* aProjectile = SelectInteractions(thePrimary);

  // now build the parton pairs. HPW
  SplitHadrons();

  // soft collisions first HPW, ordering is vital
  PerformSoftCollisions();

  // the rest is diffractive HPW
  PerformDiffractiveCollisions();

  // clean-up, if necessary
  std::for_each(theInteractions.begin(), theInteractions.end(), DeleteInteractionContent());
  theInteractions.clear();
  std::for_each(theTargets.begin(), theTargets.end(), DeleteSplitableHadron());
  theTargets.clear();
  delete aProjectile;
}

G4VSplitableHadron* G4QGSParticipants::SelectInteractions(const G4ReactionProduct  &thePrimary)
{
  G4VSplitableHadron* aProjectile = new G4QGSMSplitableHadron(thePrimary, TRUE); // @@@ check the TRUE
  G4PomeronCrossSection theProbability(thePrimary.GetDefinition()); // @@@ should be data member
  G4double outerRadius = theNucleus->GetOuterRadius();
  // Check reaction threshold

  theNucleus->StartLoop();
  G4Nucleon * pNucleon = theNucleus->GetNextNucleon();
  G4LorentzVector aPrimaryMomentum(thePrimary.GetMomentum(), thePrimary.GetTotalEnergy());
  //--DEBUG--G4cout << " qgspart- " << aPrimaryMomentum << " # " << aPrimaryMomentum.mag()
  //--DEBUG--      << pNucleon->Get4Momentum() << " # " << (pNucleon->Get4Momentum()).mag()<< G4endl;
  G4double s_nucleus = (aPrimaryMomentum + pNucleon->Get4Momentum()).mag2();
  G4double ThresholdMass = thePrimary.GetMass() + pNucleon->GetDefinition()->GetPDGMass();
  ModelMode = SOFT;
  if (sqr(ThresholdMass + ThresholdParameter) > s_nucleus)
  {
    ModelMode = DIFFRACTIVE;
    //throw G4HadronicException(__FILE__, __LINE__, 
    //                          "Initial energy is too low. The 4-vectors of the input are inconsistant with the particle masses.");
  }
  if (sqr(ThresholdMass + QGSMThreshold) > s_nucleus) // thus only diffractive in cascade!
  {
    ModelMode = DIFFRACTIVE;
  }

  // first find the collisions HPW
  std::for_each(theInteractions.begin(), theInteractions.end(), DeleteInteractionContent());
  theInteractions.clear();
  G4int totalCuts = 0;

  #ifdef debug_QGS
  G4double eK = thePrimary.GetKineticEnergy()/GeV;
  #endif
  #ifdef debug_G4QGSParticipants
  G4double impactUsed = 0;
  G4LorentzVector intNuclMom;
  #endif

  const G4int maxNumberOfLoops = 1000;
  G4int loopCounter = -1;
  while ( (theInteractions.size() == 0) && ++loopCounter < maxNumberOfLoops )  /* Loop checking, 26.10.2015, A.Ribon */
  {
    // choose random impact parameter HPW
    std::pair<G4double, G4double> theImpactParameter;
    theImpactParameter = theNucleus->ChooseImpactXandY(outerRadius+theNucleonRadius);
    G4double impactX = theImpactParameter.first;
    G4double impactY = theImpactParameter.second;

    // loop over nucleons to find collisions
    theNucleus->StartLoop();
    G4int nucleonCount = 0; // debug
    G4QGSParticipants_NPart = 0;
    #ifdef debug_G4QGSParticipants
    intNuclMom=aPrimaryMomentum;
    #endif
    while( (pNucleon = theNucleus->GetNextNucleon()) )  /* Loop checking, 26.10.2015, A.Ribon */
    {
      if(totalCuts>1.5*thePrimary.GetKineticEnergy()/GeV)
      {
        break;
      }
      nucleonCount++; // debug
      // Needs to be moved to Probability class @@@
      G4LorentzVector nucleonMomentum=pNucleon->Get4Momentum();
      nucleonMomentum.setE(nucleonMomentum.e()-pNucleon->GetBindingEnergy());
      G4double s_nucleon = (aPrimaryMomentum + nucleonMomentum).mag2();
      G4double Distance2 = sqr(impactX - pNucleon->GetPosition().x()) +
			   sqr(impactY - pNucleon->GetPosition().y());
      G4double Probability = theProbability.GetInelasticProbability(s_nucleon, Distance2);
      // test for inelastic collision
      G4double rndNumber = G4UniformRand();
      //      ModelMode = DIFFRACTIVE;
      if (Probability > rndNumber)
      {
        #ifdef debug_G4QGSParticipants
        G4cout << "DEBUG p="<< Probability<<" r="<<rndNumber<<" d="<<std::sqrt(Distance2)<<G4endl;
        G4cout << " qgspart+  " << aPrimaryMomentum << " # " << aPrimaryMomentum.mag()
	       << pNucleon->Get4Momentum() << " # " << (pNucleon->Get4Momentum()).mag()<< G4endl;
        intNuclMom += nucleonMomentum;
        #endif
        pNucleon->SetMomentum(nucleonMomentum);
	G4QGSMSplitableHadron* aTarget = new G4QGSMSplitableHadron(*pNucleon);
	G4QGSParticipants_NPart++;
	theTargets.push_back(aTarget);
	pNucleon->Hit(aTarget);
	if ((theProbability.GetDiffractiveProbability(s_nucleon, Distance2)/Probability > G4UniformRand()
	    &&(ModelMode==SOFT)) || (ModelMode==DIFFRACTIVE ))
	{
	  // diffractive interaction occurs
	  if(IsSingleDiffractive())
	  {
	    theSingleDiffExcitation.ExciteParticipants(aProjectile, aTarget);
	  } else {
	    theDiffExcitaton.ExciteParticipants(aProjectile, aTarget);
          }
	  G4InteractionContent * aInteraction = new G4InteractionContent(aProjectile);
	  aInteraction->SetTarget(aTarget);
	  theInteractions.push_back(aInteraction);
	  aInteraction->SetNumberOfDiffractiveCollisions(1);
	  totalCuts += 1;
	} else {
	  // nondiffractive soft interaction occurs
	  // sample nCut+1 (cut Pomerons) pairs of strings can be produced
	  G4int nCut;
	  G4double * running = new G4double[nCutMax];
	  running[0] = 0;
	  for(nCut = 0; nCut < nCutMax; nCut++)
	  {
            running[nCut] = theProbability.GetCutPomeronProbability(s_nucleon, Distance2, nCut + 1);
	    if(nCut!=0) running[nCut] += running[nCut-1];
	  }
	  G4double random = running[nCutMax-1]*G4UniformRand();
	  for(nCut = 0; nCut < nCutMax; nCut++)
          {
            if(running[nCut] > random) break;
	  }
	  delete [] running;
	  nCut = 0;
	  aTarget->IncrementCollisionCount(nCut+1);
	  aProjectile->IncrementCollisionCount(nCut+1);
	  G4InteractionContent * aInteraction = new G4InteractionContent(aProjectile);
	  aInteraction->SetTarget(aTarget);
	  aInteraction->SetNumberOfSoftCollisions(nCut+1);
          theInteractions.push_back(aInteraction);
	  totalCuts += nCut+1;
	  #ifdef debug_G4QGSParticipants
	  impactUsed=Distance2;
	  #endif
	}
      }
    }

    #ifdef debug_G4QGSParticipants
    G4cout << G4endl<<"NUCLEONCOUNT "<<nucleonCount<<G4endl;
    G4cout << " Interact 4-Vect " << intNuclMom << G4endl;
    #endif

  }

  if ( loopCounter >= maxNumberOfLoops ) {
    G4ExceptionDescription ed;
    ed << " loopCounter exceeds maxNumberOfLoops : forced exit! " << G4endl;
    G4Exception( "G4QGSParticipants::SelectInteractions ", "HAD_QGS_001", JustWarning, ed );
  }

  #ifdef debug_G4QGSParticipants
  G4cout << G4endl<<"CUTDEBUG "<< totalCuts <<G4endl;
  G4cout << "Impact Parameter used = "<<impactUsed<<G4endl;
  #endif

  return aProjectile;
}

void G4QGSParticipants::PerformDiffractiveCollisions()
{
  // remove the "G4PartonPair::PROJECTILE", etc., which are not necessary. @@@
  unsigned int i;
  for(i = 0; i < theInteractions.size(); i++)
  {
    G4InteractionContent* anIniteraction = theInteractions[i];
    G4VSplitableHadron* aProjectile = anIniteraction->GetProjectile();
    G4Parton* aParton = aProjectile->GetNextParton();
    G4PartonPair * aPartonPair;
    // projectile first HPW
    if (aParton)
    {
      aPartonPair = new G4PartonPair(aParton, aProjectile->GetNextAntiParton(),
				     G4PartonPair::DIFFRACTIVE, G4PartonPair::PROJECTILE);

      #ifdef debug_G4QGSPart_PDiffColl
      G4cout << "DiffPair Pro " << aPartonPair->GetParton1()->GetPDGcode() << " "
				<< aPartonPair->GetParton1()->Get4Momentum() << " "
				<< aPartonPair->GetParton1()->GetX() << " " << G4endl;
      G4cout << "         " << aPartonPair->GetParton2()->GetPDGcode() << " "
			    << aPartonPair->GetParton2()->Get4Momentum() << " "
			    << aPartonPair->GetParton2()->GetX() << " " << G4endl;
      #endif

      thePartonPairs.push_back(aPartonPair);
    }
    // then target HPW
    G4VSplitableHadron* aTarget = anIniteraction->GetTarget();
    aParton = aTarget->GetNextParton();
    if (aParton)
    {
      aPartonPair = new G4PartonPair(aParton, aTarget->GetNextAntiParton(),
				     G4PartonPair::DIFFRACTIVE, G4PartonPair::TARGET);

      #ifdef debug_G4QGSPart_PDiffColl
      G4cout << "DiffPair Tgt " << aPartonPair->GetParton1()->GetPDGcode() << " "
			        << aPartonPair->GetParton1()->Get4Momentum() << " "
				<< aPartonPair->GetParton1()->GetX() << " " << G4endl;
      G4cout << "         " << aPartonPair->GetParton2()->GetPDGcode() << " "
			    << aPartonPair->GetParton2()->Get4Momentum() << " "
			    << aPartonPair->GetParton2()->GetX() << " " << G4endl;
      #endif

      thePartonPairs.push_back(aPartonPair);
    }
  }
}

void G4QGSParticipants::PerformSoftCollisions()
{
  std::vector<G4InteractionContent*>::iterator i;
  G4LorentzVector str4Mom;
  i = theInteractions.begin();
  while ( i != theInteractions.end() )  /* Loop checking, 10.08.2015, A.Ribon */
  {
    G4InteractionContent* anIniteraction = *i;
    G4PartonPair * aPair = NULL;
    if (anIniteraction->GetNumberOfSoftCollisions())
    {
      G4VSplitableHadron* pProjectile = anIniteraction->GetProjectile();
      G4VSplitableHadron* pTarget     = anIniteraction->GetTarget();
      for (G4int j = 0; j < anIniteraction->GetNumberOfSoftCollisions(); j++)
      {
	aPair = new G4PartonPair(pTarget->GetNextParton(), pProjectile->GetNextAntiParton(),
				 G4PartonPair::SOFT, G4PartonPair::TARGET);

	#ifdef debug_G4QGSPart_PSoftColl
	G4cout << "SoftPair " << aPair->GetParton1()->GetPDGcode() << " "
			      << aPair->GetParton1()->Get4Momentum() << " "
			      << aPair->GetParton1()->GetX() << " " << G4endl;
	G4cout << "         " << aPair->GetParton2()->GetPDGcode() << " "
			      << aPair->GetParton2()->Get4Momentum() << " "
			      << aPair->GetParton2()->GetX() << " " << G4endl;
	#endif
	#ifdef debug_G4QGSParticipants
	str4Mom += aPair->GetParton1()->Get4Momentum();
	str4Mom += aPair->GetParton2()->Get4Momentum();
	#endif

        thePartonPairs.push_back(aPair);
	aPair = new G4PartonPair(pProjectile->GetNextParton(), pTarget->GetNextAntiParton(),
				 G4PartonPair::SOFT, G4PartonPair::PROJECTILE);

	#ifdef debug_G4QGSPart_PSoftColl
	G4cout << "SoftPair " << aPair->GetParton1()->GetPDGcode() << " "
	       		      << aPair->GetParton1()->Get4Momentum() << " "
			      << aPair->GetParton1()->GetX() << " " << G4endl;
        G4cout << "         " << aPair->GetParton2()->GetPDGcode() << " "
			      << aPair->GetParton2()->Get4Momentum() << " "
			      << aPair->GetParton2()->GetX() << " " << G4endl;
	#endif
	#ifdef debug_G4QGSParticipants
	str4Mom += aPair->GetParton1()->Get4Momentum();
	str4Mom += aPair->GetParton2()->Get4Momentum();
	#endif

	thePartonPairs.push_back(aPair);
      }
      delete *i;
      i=theInteractions.erase(i);    // i now points to the next interaction
    } else {
      i++;
    }
  }
 
 #ifdef debug_G4QGSPart_PSoftColl
 G4cout << " string 4 mom " << str4Mom << G4endl;
 #endif
}

