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
#ifndef G4HeavyIonParticipants_h
#define G4HeavyIonParticipants_h 1

#include "Randomize.hh"
#include "G4VParticipants.hh"
#include "G4Nucleon.hh"
#include "G4InteractionContent.hh"
#include "G4PomeronCrossSection.hh"
#include "G4DiffractiveExcitation.hh"
#include "G4SingleDiffractiveExcitation.hh"
#include "G4PartonPair.hh" 
#include "G4QGSMSplitableHadron.hh" 

G4int G4IonParticipants_NPart = 0;
G4int G4IonParticipants_NPro = 0;

class G4HeavyIonParticipants : public G4VParticipants
{
    enum  { SOFT, HARD, DIFFRACTIVE };
  public:
    G4HeavyIonParticipants();
    G4HeavyIonParticipants(const G4HeavyIonParticipants &right);
    const G4HeavyIonParticipants & operator=(const G4HeavyIonParticipants &right);
    ~G4HeavyIonParticipants(); 

    int operator==(const G4HeavyIonParticipants &right) const;
    int operator!=(const G4HeavyIonParticipants &right) const;

    G4PartonPair* GetNextPartonPair();
    void BuildInteractions(const G4ReactionProduct  &thePrimary);
    void StartPartonPairLoop();
    virtual void DoLorentzBoost(Hep3Vector aBoost) 
    {
      if(theNucleus) theNucleus->DoLorentzBoost(aBoost);
      theBoost = aBoost;
    }
    
  private:
    void SplitHadrons(); 
    void PerformSoftCollisions();
    void PerformDiffractiveCollisions();
      
  private:
    std::vector<G4InteractionContent*> theInteractions;
    struct DeleteInteractionContent {void operator()(G4InteractionContent*aC){delete aC;}};
    std::vector<G4VSplitableHadron*>   theTargets; 
    std::vector<G4VSplitableHadron*>   theProjectiles; 
    struct DeleteSplitableHadron{void operator()(G4VSplitableHadron*aS){delete aS;}};
    std::vector<G4PartonPair*>   thePartonPairs;
    struct DeletePartonPair{void operator()(G4PartonPair*aP){delete aP;}};
    G4Fancy3DNucleus theIncoming;
  
    G4SingleDiffractiveExcitation theSingleDiffExcitation;
    G4DiffractiveExcitation theDiffExcitaton;
    G4int ModelMode;
    G4bool IsSingleDiffractive();
    
    G4ThreeVector theBoost;
 
  private:
    // model parameters HPW
    const G4int nCutMax; 
    const G4double ThersholdParameter; 
    const G4double QGSMThershold; 
    const G4double theNucleonRadius;
    G4PomeronCrossSection theProtonProbability;
    G4PomeronCrossSection theNeutronProbability;
};


inline G4bool G4HeavyIonParticipants::IsSingleDiffractive()
{
  G4bool result=false;
  if(G4UniformRand()<1.) result = true;
  return result;
}

inline void G4HeavyIonParticipants::StartPartonPairLoop()
{
}

inline G4PartonPair* G4HeavyIonParticipants::GetNextPartonPair()
{
  if (thePartonPairs.empty()) return 0;
  G4PartonPair * result = thePartonPairs.back();
  thePartonPairs.pop_back();
  return result;
}


inline void G4HeavyIonParticipants::SplitHadrons()
{
  unsigned int i;
  for(i = 0; i < theInteractions.size(); i++) 
  {
    theInteractions[i]->SplitHadrons();
  }
}

#endif


