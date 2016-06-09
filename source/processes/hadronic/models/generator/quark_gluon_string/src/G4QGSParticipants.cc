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
#include "G4QGSParticipants.hh"
#include "G4LorentzVector.hh"
#include "G4Pair.hh"

// Class G4QGSParticipants 

// HPW Feb 1999
// Promoting model parameters from local variables class properties

G4QGSParticipants::G4QGSParticipants() : theDiffExcitaton(0.7*GeV, 250*MeV, 250*MeV),
                                         nCutMax(7),ThersholdParameter(0.45*GeV),
                                         QGSMThershold(3*GeV),theNucleonRadius(1.5*fermi)
                                         
{
}

G4QGSParticipants::G4QGSParticipants(const G4QGSParticipants &right)
: G4VParticipants(), nCutMax(right.nCutMax),ThersholdParameter(right.ThersholdParameter),
  QGSMThershold(right.QGSMThershold),theNucleonRadius(right.theNucleonRadius)
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
  G4double s = (aPrimaryMomentum + pNucleon->Get4Momentum()).mag2();
  G4double ThresholdMass = thePrimary.GetMass() + pNucleon->GetDefinition()->GetPDGMass(); 
  ModelMode = SOFT;
  if (sqr(ThresholdMass + ThersholdParameter) > s)
  {
    G4Exception("Initial energy is too low. The 4-vectors of the input are inconsistant with the particle masses.");
  }
  if (sqr(ThresholdMass + QGSMThershold) > s) // thus only diffractive in cascade!
  {
    ModelMode = DIFFRACTIVE;
  }
 
  // first find the collisions HPW
  std::for_each(theInteractions.begin(), theInteractions.end(), DeleteInteractionContent());
  theInteractions.clear();
  G4int totalCuts = 0;
  G4double impactUsed = 0;

  #ifdef debug_QGS
  G4double eK = thePrimary.GetKineticEnergy()/GeV;
  #endif

  while(theInteractions.size() == 0)
  {
    // choose random impact parameter HPW
    G4Pair<G4double, G4double> theImpactParameter;
    theImpactParameter = theNucleus->ChooseImpactXandY(outerRadius+theNucleonRadius);
    G4double impactX = theImpactParameter.first; 
    G4double impactY = theImpactParameter.second;
    
    // loop over nuclei to find collissions HPW
    theNucleus->StartLoop();
    G4int nucleonCount = 0; // debug
    while( (pNucleon = theNucleus->GetNextNucleon()) )
    {
      if(totalCuts>1.5*thePrimary.GetKineticEnergy()/GeV) 
      {
	break;
      }
       nucleonCount++; // debug
      // Needs to be moved to Probability class @@@
      G4double s = (aPrimaryMomentum + pNucleon->Get4Momentum()).mag2();
      G4double Distance2 = sqr(impactX - pNucleon->GetPosition().x()) +
                           sqr(impactY - pNucleon->GetPosition().y());
      G4double Probability = theProbability.GetInelasticProbability(s, Distance2);  
      // test for inelastic collision
      G4double rndNumber = G4UniformRand();
//      ModelMode = DIFFRACTIVE;
      if (Probability > rndNumber)
      {
//--DEBUG--        cout << "DEBUG p="<< Probability<<" r="<<rndNumber<<" d="<<sqrt(Distance2)<<G4endl;
        G4QGSMSplitableHadron* aTarget = new G4QGSMSplitableHadron(*pNucleon);
        theTargets.push_back(aTarget);
 	pNucleon->Hit(aTarget);
        if ((theProbability.GetDiffractiveProbability(s, Distance2)/Probability > G4UniformRand() 
             &&(ModelMode==SOFT)) || (ModelMode==DIFFRACTIVE ))
	{ 
	  // diffractive interaction occurs
	  if(IsSingleDiffractive())
	  {
	    theSingleDiffExcitation.ExciteParticipants(aProjectile, aTarget);
	  }
	  else
	  {
	    theDiffExcitaton.ExciteParticipants(aProjectile, aTarget);
	  }
          G4InteractionContent * aInteraction = new G4InteractionContent(aProjectile);
          aInteraction->SetTarget(aTarget); 
          theInteractions.push_back(aInteraction);
	  aInteraction->SetNumberOfDiffractiveCollisions(1);
          totalCuts += 1;
	}
	else
	{
	  // nondiffractive soft interaction occurs
	  // sample nCut+1 (cut Pomerons) pairs of strings can be produced
          G4int nCut;
          G4double * running = new G4double[nCutMax];
          running[0] = 0;
          for(nCut = 0; nCut < nCutMax; nCut++)
          {
	    running[nCut] = theProbability.GetCutPomeronProbability(s, Distance2, nCut + 1);
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
          impactUsed=Distance2;
       }
      }
    }
//--DEBUG--  cout << G4endl<<"NUCLEONCOUNT "<<nucleonCount<<G4endl;
  }
//--DEBUG--  cout << G4endl<<"CUTDEBUG "<< totalCuts <<G4endl;
//--DEBUG--  cout << "Impact Parameter used = "<<impactUsed<<G4endl;
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
                                     G4PartonPair::DIFFRACTIVE, 
                                     G4PartonPair::PROJECTILE);
      thePartonPairs.push_back(aPartonPair);
    }
    // then target HPW
    G4VSplitableHadron* aTarget = anIniteraction->GetTarget();
    aParton = aTarget->GetNextParton();
    if (aParton)
    {
      aPartonPair = new G4PartonPair(aParton, aTarget->GetNextAntiParton(), 
                                     G4PartonPair::DIFFRACTIVE, 
                                     G4PartonPair::TARGET);
      thePartonPairs.push_back(aPartonPair);
    }
  }
}

void G4QGSParticipants::PerformSoftCollisions()
{
  std::vector<G4InteractionContent*>::iterator i;
  for(i = theInteractions.begin(); i != theInteractions.end(); i++)   
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
        thePartonPairs.push_back(aPair);
        aPair = new G4PartonPair(pProjectile->GetNextParton(), pTarget->GetNextAntiParton(), 
                                 G4PartonPair::SOFT, G4PartonPair::PROJECTILE);
        thePartonPairs.push_back(aPair);
      }  
      delete *i;
      i=theInteractions.erase(i);
      i--;
    }
  }
}
