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
#include "G4HeavyIonParticipants.hh"
#include "G4LorentzVector.hh"
#include "G4Pair.hh"
#include "g4std/algorithm"

// Class G4HeavyIonParticipants 

// HPW April 2001
// Promoting model parameters from local variables class properties

G4HeavyIonParticipants::G4HeavyIonParticipants() : theDiffExcitaton(0.7*GeV, 250*MeV, 250*MeV),
                                         nCutMax(7),ThersholdParameter(0.45*GeV),
                                         QGSMThershold(3*GeV),theNucleonRadius(1.5*fermi)
{
}

G4HeavyIonParticipants::G4HeavyIonParticipants(const G4HeavyIonParticipants &right)
: nCutMax(right.nCutMax),ThersholdParameter(right.ThersholdParameter),
  QGSMThershold(right.QGSMThershold),theNucleonRadius(right.theNucleonRadius)
{
}


G4HeavyIonParticipants::~G4HeavyIonParticipants()
{
}

void G4HeavyIonParticipants::BuildInteractions(const G4ReactionProduct  &thePrimary) 
{
  G4PomeronCrossSection theProtonProbability(G4Proton::Proton()); // @@@ should be data member
  G4PomeronCrossSection theNeutronProbability(G4Neutron::Neutron()); // @@@ should be data member
  G4double outerRadiusA = theNucleus->GetOuterRadius();
  theIncoming.Init(thePrimary.GetDefinition()->GetBaryonNumber(),
                   thePrimary.GetDefinition()->GetPDGCharge() );
  theIncoming.DoLorentzBoost(-theBoost);
  G4double outerRadiusB = theIncoming.GetOuterRadius();
  ModelMode = SOFT;
 
  // first find the collisions HPW
  G4std::for_each(theInteractions.begin(), theInteractions.end(), DeleteInteractionContent());
  theInteractions.clear();
  G4int totalCuts = 0;
  G4double impactUsed = 0;
  G4int deb1, deb2;
  while(theInteractions.empty())
  {
    // choose random impact parameter HPW
    G4Pair<G4double, G4double> theImpactParameter;
    theImpactParameter = theNucleus->ChooseImpactXandY(outerRadiusA+outerRadiusB+theNucleonRadius);
    G4double impactX = theImpactParameter.first; 
    G4double impactY = theImpactParameter.second;
    
    // loop over nuclei to find collissions HPW
    theIncoming.StartLoop();
    G4Nucleon * iNucleon;
    G4Nucleon * pNucleon;
    deb1 = 0;
    while( (iNucleon = theIncoming.GetNextNucleon()) )
    {
      deb1 ++;
      deb2 = 0;
      theNucleus->StartLoop();
      G4int nucleonCount = 0; // debug
      while( (pNucleon = theNucleus->GetNextNucleon()) )
      {
        deb2 ++;
	if(totalCuts>nCutMax) break;
         nucleonCount++; // debug
        // Needs to be moved to Probability class @@@
        G4double s = (iNucleon->Get4Momentum() + pNucleon->Get4Momentum()).mag2();
        G4double Distance2 = sqr(impactX - pNucleon->GetPosition().x() - iNucleon->GetPosition().x()) +
                             sqr(impactY - pNucleon->GetPosition().y() - iNucleon->GetPosition().y());
        G4PomeronCrossSection * theProbability=NULL;
        G4double Probability;
        if(iNucleon->GetDefinition() == G4Proton::Proton()) theProbability = &theProtonProbability;
        if(iNucleon->GetDefinition() == G4Neutron::Neutron()) theProbability = &theNeutronProbability;  
        Probability = theNeutronProbability.GetInelasticProbability(s, Distance2);
// test for inelastic collision
        G4double rndNumber = G4UniformRand();
//      ModelMode = DIFFRACTIVE;
//        if(deb1<2) 
//	{
//	  Probability = 5; // HPW for debugging
//        }
//	else
//	{
//	  Probability = -5;
//	}
	if (Probability > rndNumber)
        {

//--DEBUG--        cout << "DEBUG p="<< Probability<<" r="<<rndNumber<<" d="<<sqrt(Distance2)<<G4endl;
          G4VSplitableHadron* aTarget;
	  if(pNucleon->GetSplitableHadron())
	  {
	    aTarget = pNucleon->GetSplitableHadron();
	  }
	  else 
	  {
	    aTarget = new G4QGSMSplitableHadron(*pNucleon);
            theTargets.push_back(aTarget);
 	    pNucleon->Hit(aTarget);
	  }
          G4VSplitableHadron* aProject;
	  if(iNucleon->GetSplitableHadron())
	  {
	    aProject = iNucleon->GetSplitableHadron();
	  }
	  else
	  {
	    aProject = new G4QGSMSplitableHadron(*iNucleon, TRUE);
	    theProjectiles.push_back(aProject);
 	    iNucleon->Hit(aProject);
          }
          if ((theProbability->GetDiffractiveProbability(s, Distance2)/Probability > G4UniformRand() 
               &&(ModelMode==SOFT)) || (ModelMode==DIFFRACTIVE ))
        	{ 
	          // diffractive interaction occurs
	          if(IsSingleDiffractive())
	          {
	            theSingleDiffExcitation.ExciteParticipants(aProject, aTarget);
	          }
	          else
	          {
	            theDiffExcitaton.ExciteParticipants(aProject, aTarget);
	          }
            G4InteractionContent * aInteraction = new G4InteractionContent(aProject);
            aInteraction->SetTarget(aTarget); 
            theInteractions.push_back(aInteraction);
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
	            running[nCut] = theProbability->GetCutPomeronProbability(s, Distance2, nCut + 1);
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
            aProject->IncrementCollisionCount(nCut+1);
            G4InteractionContent * aInteraction = new G4InteractionContent(aProject);
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
  }
//--DEBUG--  cout << G4endl<<"CUTDEBUG "<< totalCuts <<G4endl;
//--DEBUG--  cout << "Impact Parameter used = "<<impactUsed<<G4endl;
  // now build the parton pairs. HPW
  SplitHadrons();
  
  // soft collisions first HPW, ordering is vital
  PerformSoftCollisions();
   
  // the rest is diffractive HPW
  PerformDiffractiveCollisions();
  
  // clean-up, if necessary
  G4std::for_each(theInteractions.begin(), theInteractions.end(), DeleteInteractionContent());
  theInteractions.clear();
  G4std::for_each(theTargets.begin(), theTargets.end(), DeleteSplitableHadron());
  theTargets.clear();
  G4std::for_each(theProjectiles.begin(), theProjectiles.end(), DeleteSplitableHadron());
  theProjectiles.clear();
}

void G4HeavyIonParticipants::PerformDiffractiveCollisions()
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

void G4HeavyIonParticipants::PerformSoftCollisions()
{
  G4std::vector<G4InteractionContent*>::iterator i;
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
