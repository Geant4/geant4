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
#include "globals.hh"
#include "G4GammaParticipants.hh"
#include "G4LorentzVector.hh"
#include "G4V3DNucleus.hh"
#include <utility>

// Class G4GammaParticipants 

// J.P. Wellisch, April 2002
// new participants class for gamma nuclear, with this design more can come with 
// cross-section based, and quasi-eiconal model based modelling 
//
// 20110805  M. Kelsey -- Follow change to G4V3DNucleus::GetNucleons()

G4VSplitableHadron* G4GammaParticipants::SelectInteractions(const G4ReactionProduct  &thePrimary) 
{
  // Check reaction threshold  - goes to CheckThreshold
  G4VSplitableHadron* aProjectile = new G4QGSMSplitableHadron(thePrimary, TRUE); // @@@ check the TRUE

  const std::vector<G4Nucleon>& theTargetNuc = theNucleus->GetNucleons();
  G4LorentzVector aPrimaryMomentum(thePrimary.GetMomentum(), thePrimary.GetTotalEnergy());
  if((!(aPrimaryMomentum.e()>-1)) && (!(aPrimaryMomentum.e()<1)) )
  {
    throw G4HadronicException(__FILE__, __LINE__,
			      "G4GammaParticipants::SelectInteractions: primary nan energy.");
  }
  G4double S = (aPrimaryMomentum + theTargetNuc[0].Get4Momentum()).mag2();
  G4double ThresholdMass = thePrimary.GetMass() + theTargetNuc[0].GetDefinition()->GetPDGMass();
  ModelMode = SOFT;
  if (sqr(ThresholdMass + ThresholdParameter) > S)
  {
    ModelMode = DIFFRACTIVE;
    //throw G4HadronicException(__FILE__, __LINE__, 
    //                          "Initial energy is too low. The 4-vectors of the input are inconsistant with the particle masses.");
  }
  if (sqr(ThresholdMass + QGSMThreshold) > S) // thus only diffractive in cascade!
  {
    ModelMode = DIFFRACTIVE;
  }

  // first find the collisions HPW
  std::for_each(theInteractions.begin(), theInteractions.end(), DeleteInteractionContent());
  theInteractions.clear();
  G4int totalCuts = 0;

  #ifdef debug_G4GammaParticipants
  G4double eK = thePrimary.GetKineticEnergy()/GeV;
  G4int nucleonCount = theTargetNuc.size(); // debug
  #endif

  G4int theCurrent = static_cast<G4int> (theTargetNuc.size()*G4UniformRand());
  const G4Nucleon& pNucleon = theTargetNuc[theCurrent];
  G4QGSMSplitableHadron* aTarget = new G4QGSMSplitableHadron(pNucleon);
  theTargets.push_back(aTarget);
  const_cast<G4Nucleon&>(pNucleon).Hit(aTarget);
  if ( (0.06 > G4UniformRand() &&(ModelMode==SOFT)) || (ModelMode==DIFFRACTIVE ) )
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

