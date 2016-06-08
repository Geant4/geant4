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
#ifndef G4QGSParticipants_h
#define G4QGSParticipants_h 1

#include "Randomize.hh"
#include "G4VParticipants.hh"
#include "G4Nucleon.hh"
#include "G4InteractionContent.hh"
#include "G4PomeronCrossSection.hh"
#include "G4DiffractiveExcitation.hh"
#include "G4SingleDiffractiveExcitation.hh"
#include "G4PartonPair.hh" 
#include "G4QGSMSplitableHadron.hh" 

class G4QGSParticipants : public G4VParticipants
{
    enum  { SOFT, HARD, DIFFRACTIVE };
  public:
    G4QGSParticipants();
    G4QGSParticipants(const G4QGSParticipants &right);
    const G4QGSParticipants & operator=(const G4QGSParticipants &right);
    ~G4QGSParticipants(); 

    int operator==(const G4QGSParticipants &right) const;
    int operator!=(const G4QGSParticipants &right) const;

    G4PartonPair* GetNextPartonPair();
    void BuildInteractions(const G4ReactionProduct  &thePrimary);
    void StartPartonPairLoop();
    
  private:
    void SplitHadrons(); 
    void PerformSoftCollisions();
    void PerformDiffractiveCollisions();
      
  private:
    G4RWTPtrOrderedVector<G4InteractionContent> theInteractions;
    G4RWTPtrOrderedVector<G4VSplitableHadron>   theTargets; 
    G4RWTPtrOrderedVector<G4PartonPair>   thePartonPairs;
  
    G4SingleDiffractiveExcitation theSingleDiffExcitation;
    G4DiffractiveExcitation theDiffExcitaton;
    G4int ModelMode;
    G4bool IsSingleDiffractive();
 
  private:
    // model parameters HPW
    const G4int nCutMax; 
    const G4double ThersholdParameter; 
    const G4double QGSMThershold; 
    const G4double theNucleonRadius;
    
};


inline G4bool G4QGSParticipants::IsSingleDiffractive()
{
  G4bool result;
  if(G4UniformRand()<1.) result = true;
  return result;
}

inline void G4QGSParticipants::StartPartonPairLoop()
{
}

inline G4PartonPair* G4QGSParticipants::GetNextPartonPair()
{
if (thePartonPairs.isEmpty()) return 0;
return thePartonPairs.removeLast();
}


inline void G4QGSParticipants::SplitHadrons()
{
  G4int i;
  for(i = 0; i < theInteractions.length(); i++) 
  {
    theInteractions.at(i)->SplitHadrons();
  }
}

#endif


