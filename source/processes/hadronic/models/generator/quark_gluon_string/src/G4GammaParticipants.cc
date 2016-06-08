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
#include "G4GammaParticipants.hh"
#include "G4LorentzVector.hh"
#include "G4Pair.hh"

// Class G4GammaParticipants 

// J.P. Wellisch, April 2002
// new participants class for gamma nuclear, with this design more can come with 
// cross-section based, and quasi-eiconal model based modelling 

G4VSplitableHadron* G4GammaParticipants::SelectInteractions(const G4ReactionProduct  &thePrimary) 
{
  // Check reaction threshold  - goes to CheckThreshold
  G4VSplitableHadron* aProjectile = new G4QGSMSplitableHadron(thePrimary, TRUE); // @@@ check the TRUE

  const G4std::vector<G4Nucleon *> & theTargetNuc = theNucleus->GetNucleons();
  G4LorentzVector aPrimaryMomentum(thePrimary.GetMomentum(), thePrimary.GetTotalEnergy());
  G4double s = (aPrimaryMomentum + theTargetNuc[0]->Get4Momentum()).mag2();
  G4double ThresholdMass = thePrimary.GetMass() + theTargetNuc[0]->GetDefinition()->GetPDGMass(); 
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
  G4std::for_each(theInteractions.begin(), theInteractions.end(), DeleteInteractionContent());
  theInteractions.clear();
  G4int totalCuts = 0;

   #ifdef debug_G4GammaParticipants
   G4double eK = thePrimary.GetKineticEnergy()/GeV;
   G4int nucleonCount = theTargetNuc.size(); // debug
   #endif

  G4int theCurrent = static_cast<G4int> (theTargetNuc.size()*G4UniformRand());
  G4Nucleon * pNucleon = theTargetNuc[theCurrent];
  G4QGSMSplitableHadron* aTarget = new G4QGSMSplitableHadron(*pNucleon);
  theTargets.push_back(aTarget);
   pNucleon->Hit(aTarget);
   if ( (0.06 > G4UniformRand() &&(ModelMode==SOFT)) || (ModelMode==DIFFRACTIVE ) )
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
     aTarget->IncrementCollisionCount(1);
     aProjectile->IncrementCollisionCount(1);
     G4InteractionContent * aInteraction = new G4InteractionContent(aProjectile);
     aInteraction->SetTarget(aTarget);
     aInteraction->SetNumberOfSoftCollisions(1);
     theInteractions.push_back(aInteraction);
     totalCuts += 1;
    }
  return aProjectile;
}
