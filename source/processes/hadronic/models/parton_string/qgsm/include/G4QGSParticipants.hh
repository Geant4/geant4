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
#ifndef G4QGSParticipants_h
#define G4QGSParticipants_h 1

#include "Randomize.hh"
#include "G4VParticipants.hh"
#include "G4Nucleon.hh"
#include "G4InteractionContent.hh"
#include "G4PomeronCrossSection.hh"
#include "G4QGSDiffractiveExcitation.hh"
#include "G4SingleDiffractiveExcitation.hh"
#include "G4PartonPair.hh" 
#include "G4QGSMSplitableHadron.hh" 
#include "G4V3DNucleus.hh"


class G4QGSParticipants : public G4VParticipants
{
  public:
    G4QGSParticipants();
    G4QGSParticipants(const G4QGSParticipants &right);
    const G4QGSParticipants & operator=(const G4QGSParticipants &right);
    virtual ~G4QGSParticipants();

    int operator==(const G4QGSParticipants &right) const;
    int operator!=(const G4QGSParticipants &right) const;

    virtual void DoLorentzBoost(G4ThreeVector aBoost)
    {
      if(theNucleus) theNucleus->DoLorentzBoost(aBoost);
      theBoost = aBoost;
    }

    G4PartonPair* GetNextPartonPair();
    void BuildInteractions(const G4ReactionProduct  &thePrimary);
    void StartPartonPairLoop();

  protected:
    virtual G4VSplitableHadron* SelectInteractions(const G4ReactionProduct  &thePrimary);
    void SplitHadrons(); 
    void PerformSoftCollisions();
    void PerformDiffractiveCollisions();

  protected:
    struct DeleteInteractionContent {void operator()(G4InteractionContent*aC){delete aC;}};
    std::vector<G4InteractionContent*> theInteractions;
    struct DeleteSplitableHadron{void operator()(G4VSplitableHadron*aS){delete aS;}};
    std::vector<G4VSplitableHadron*>   theTargets;
    struct DeletePartonPair{void operator()(G4PartonPair*aP){delete aP;}};
    std::vector<G4PartonPair*>   thePartonPairs;

    G4SingleDiffractiveExcitation theSingleDiffExcitation;
    G4QGSDiffractiveExcitation theDiffExcitaton;
    G4int ModelMode;
    G4bool IsSingleDiffractive();

    G4ThreeVector theBoost;

  protected:
    // model parameters HPW
    enum  { SOFT, DIFFRACTIVE };
    const G4int nCutMax;
    const G4double ThresholdParameter;
    const G4double QGSMThreshold;
    const G4double theNucleonRadius;
};


inline G4bool G4QGSParticipants::IsSingleDiffractive()
{
  G4bool result=false;
  if(G4UniformRand()<1.) result = true;
  return result;
}

inline void G4QGSParticipants::StartPartonPairLoop()
{
}

inline G4PartonPair* G4QGSParticipants::GetNextPartonPair()
{
  if (thePartonPairs.empty()) return 0;
  G4PartonPair * result = thePartonPairs.back();
  thePartonPairs.pop_back();
  return result;
}

inline void G4QGSParticipants::SplitHadrons()
{
  unsigned int i;
  for(i = 0; i < theInteractions.size(); i++)
  {
    theInteractions[i]->SplitHadrons();
  }
}

#endif

