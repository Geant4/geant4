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
//  Class G4QMDFermiParticipants 
//  Participants for nucleus projectile  
//  Created by T. Koi SLAC/SCCS 
//  Strongly based on G4QGSParticipants
//
#include "globals.hh"
#include "G4QMDFermiParticipants.hh"
#include "G4LorentzVector.hh"
#include <utility>

// HPW Feb 1999
// Promoting model parameters from local variables class properties
G4int G4QMDFermiParticipants_NPart = 0;

G4QMDFermiParticipants::G4QMDFermiParticipants() :G4QGSParticipants() ,theDiffExcitaton(), //0.7*GeV, 250*MeV, 250*MeV),
                                         nCutMax(7),ThersholdParameter(0.45*GeV),
                                         QGSMThershold(3*GeV),theNucleonRadius(1.5*fermi)
                                         
{
}

G4QMDFermiParticipants::G4QMDFermiParticipants(const G4QMDFermiParticipants &right)
: G4QGSParticipants(), nCutMax(right.nCutMax),ThersholdParameter(right.ThersholdParameter),
  QGSMThershold(right.QGSMThershold),theNucleonRadius(right.theNucleonRadius)
{
}


G4QMDFermiParticipants::~G4QMDFermiParticipants()
{
}

void G4QMDFermiParticipants::BuildInteractions(const G4ReactionProduct  &thePrimary) 
{
  
  // Find the collisions and collition conditions
//081116
//G4cout << "G4QMDFermiParticipants::BuildInteractions " << G4endl;
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
//TK 
  std::for_each(theProjs.begin(), theProjs.end(), DeleteSplitableHadron());
  theProjs.clear();
//TK 
  delete aProjectile;

}

G4VSplitableHadron* G4QMDFermiParticipants::SelectInteractions(const G4ReactionProduct  &thePrimary)
{
  //TK
  G4VSplitableHadron* axProjectile = new G4QGSMSplitableHadron(thePrimary, TRUE); // @@@ check the TRUE
  //G4PomeronCrossSection theProbability(thePrimary.GetDefinition()); // @@@ should be data member
  G4double outerRadius = theNucleus->GetOuterRadius();
  // Check reaction threshold 

  theNucleus->StartLoop();
  G4Nucleon * pNucleon = theNucleus->GetNextNucleon();
  G4LorentzVector aPrimaryMomentum(thePrimary.GetMomentum(), thePrimary.GetTotalEnergy());
//--DEBUG--G4cout << " qgspart- " << aPrimaryMomentum << " # " << aPrimaryMomentum.mag() 
//--DEBUG--      << pNucleon->Get4Momentum() << " # " << (pNucleon->Get4Momentum()).mag()<< G4endl;
  G4double s = (aPrimaryMomentum + pNucleon->Get4Momentum()).mag2();
  G4double ThresholdMass = thePrimary.GetMass() + pNucleon->GetDefinition()->GetPDGMass(); 
  ModelMode = SOFT;
  if (sqr(ThresholdMass + ThersholdParameter) > s)
  {
    throw G4HadronicException(__FILE__, __LINE__, "Initial energy is too low. The 4-vectors of the input are inconsistant with the particle masses.");
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
#ifdef debug_G4QMDFermiParticipants
  G4LorentzVector intNuclMom;
#endif
//TK
   G4Fancy3DNucleus* proj = NULL;  
   G4Nucleon* projNucleon = NULL; 
   proj = new G4Fancy3DNucleus();  
   proj->Init( thePrimary.GetDefinition()->GetAtomicMass() , thePrimary.GetDefinition()->GetAtomicNumber() ); 
   proj->SortNucleonsInZ();
   proj->StartLoop();
   proj->CenterNucleons();
/*
   while ( ( projNucleon = proj->GetNextNucleon() ) )
   {
      G4cout << projNucleon->GetParticleType()->GetParticleName() << " " << projNucleon->GetPosition()/fermi << G4endl;
   } 
   G4cout << proj->GetOuterRadius()/fermi << G4endl;
*/

//TK
   while( theInteractions.size() == 0 )
   {
       // choose random impact parameter HPW
       std::pair<G4double, G4double> theImpactParameter;
       theImpactParameter = theNucleus->ChooseImpactXandY(outerRadius+ proj->GetOuterRadius() );

       G4double impactX = theImpactParameter.first; 
       G4double impactY = theImpactParameter.second;

       //G4cout << "Impact "  << impactX/fermi << " " << impactY/fermi << G4endl;
    
    G4int nucleonCount = 0; // debug
    G4QMDFermiParticipants_NPart = 0;
//TK
    proj->StartLoop();
    while ( projNucleon = proj->GetNextNucleon() )
    {

       impactX = theImpactParameter.first + projNucleon->GetPosition().x();
       impactY = theImpactParameter.second + projNucleon->GetPosition().y();


        // for p or n result is same 
        G4VSplitableHadron* aProjectile = new G4QGSMSplitableHadron(thePrimary, TRUE); // @@@ check the TRUE
        theProjs.push_back( aProjectile );
        G4PomeronCrossSection theProbability( G4Proton::Proton() ); 

        aProjectile->SetDefinition( projNucleon->GetParticleType() ); 
        aProjectile->Set4Momentum( aPrimaryMomentum/thePrimary.GetDefinition()->GetAtomicMass() ); 

    // loop over nuclei to find collissions HPW
    theNucleus->StartLoop();
#ifdef debug_G4QMDFermiParticipants
    intNuclMom=aPrimaryMomentum;
#endif
    while( (pNucleon = theNucleus->GetNextNucleon()) )
    {
if ( pNucleon->AreYouHit() ) break;
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
//--DEBUG--        cout << "DEBUG p="<< Probability<<" r="<<rndNumber<<" d="<<std::sqrt(Distance2)<<G4endl;
//--DEBUG--         G4cout << " qgspart+  " << aPrimaryMomentum << " # " << aPrimaryMomentum.mag() 
//--DEBUG--              << pNucleon->Get4Momentum() << " # " << (pNucleon->Get4Momentum()).mag()<< G4endl;

#ifdef debug_G4QMDFermiParticipants
       intNuclMom += pNucleon->Get4Momentum();
#endif
       G4QGSMSplitableHadron* aTarget = new G4QGSMSplitableHadron(*pNucleon);
        G4QMDFermiParticipants_NPart ++;
	theTargets.push_back(aTarget);
 	pNucleon->Hit(aTarget);

        //TK
        projNucleon->Hit(aTarget);
        //G4cout << "Hit " << projNucleon->GetParticleType()->GetParticleName() << G4endl;
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

    } // targ loop
//--DEBUG--  g4cout << G4endl<<"NUCLEONCOUNT "<<nucleonCount<<G4endl;
//--DEBUG--  G4cout << " Interact 4-Vect " << intNuclMom << G4endl; 
//TK
//G4cout << "projNucleon->AreYouHit() " << projNucleon->AreYouHit() << G4endl;
  } // proj loop
  }
   //TK
   delete proj;
//--DEBUG--  cout << G4endl<<"CUTDEBUG "<< totalCuts <<G4endl;
//--DEBUG--  cout << "Impact Parameter used = "<<impactUsed<<G4endl;
  //TK
  return axProjectile;
}



void G4QMDFermiParticipants::PerformDiffractiveCollisions()
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
#ifdef debug_G4QMDFermiPart_PDiffColl
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
                                     G4PartonPair::DIFFRACTIVE, 
                                     G4PartonPair::TARGET);
#ifdef debug_G4QMDFermiPart_PDiffColl
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


void G4QMDFermiParticipants::PerformSoftCollisions()
{
  std::vector<G4InteractionContent*>::iterator i;
  G4LorentzVector str4Mom;
  i = theInteractions.begin();
  while ( i != theInteractions.end() )   
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
#ifdef debug_G4QMDFermiPart_PSoftColl
	G4cout << "SoftPair " << aPair->GetParton1()->GetPDGcode() << " " 
			      << aPair->GetParton1()->Get4Momentum() << " "
			      << aPair->GetParton1()->GetX() << " " << G4endl;
	G4cout << "         " << aPair->GetParton2()->GetPDGcode() << " " 
			      << aPair->GetParton2()->Get4Momentum() << " "
			      << aPair->GetParton2()->GetX() << " " << G4endl;
#endif
#ifdef debug_G4QMDFermiParticipants
        str4Mom += aPair->GetParton1()->Get4Momentum();
	str4Mom += aPair->GetParton2()->Get4Momentum();
#endif
        thePartonPairs.push_back(aPair);
        aPair = new G4PartonPair(pProjectile->GetNextParton(), pTarget->GetNextAntiParton(), 
                                 G4PartonPair::SOFT, G4PartonPair::PROJECTILE);
#ifdef debug_G4QMDFermiPart_PSoftColl
	G4cout << "SoftPair " << aPair->GetParton1()->GetPDGcode() << " " 
			      << aPair->GetParton1()->Get4Momentum() << " "
			      << aPair->GetParton1()->GetX() << " " << G4endl;
	G4cout << "         " << aPair->GetParton2()->GetPDGcode() << " " 
			      << aPair->GetParton2()->Get4Momentum() << " "
			      << aPair->GetParton2()->GetX() << " " << G4endl;
#endif
#ifdef debug_G4QMDFermiParticipants
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
#ifdef debug_G4QMDFermiPart_PSoftColl
  G4cout << " string 4 mom " << str4Mom << G4endl;
#endif
}
