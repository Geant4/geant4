// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleChangeForRadDecay.hh,v 1.2 2000-02-13 15:08:20 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
// 
// ------------------------------------------------------------
//   Implemented for the new scheme                 25 Jan. 2000  H.Kurahige
//
// Class Description
//  This class is a concrete class for ParticleChange for RadDecay
//        
#ifndef G4ParticleChangeForRadDecay_h
#define G4ParticleChangeForRadDecay_h 1

#include "globals.hh"
#include "G4ios.hh"
class G4VTouchable;
#include "G4ParticleChange.hh"

class G4ParticleChangeForRadDecay: public G4ParticleChange
{ 
  public:
    // default constructor
    G4ParticleChangeForRadDecay(){}

    // destructor
    virtual ~G4ParticleChangeForRadDecay(){}

  protected:
    // hide copy constructor and assignment operaor as protected
    G4ParticleChangeForRadDecay(const G4ParticleChangeForRadDecay &right){}
    G4ParticleChangeForRadDecay & operator=(const G4ParticleChangeForRadDecay &right){}

  public: // with description
  void AddSecondary(G4Track* aSecondary);
    //  Add a secondary particle to theListOfSecondaries.
    // ------------------------------------------------------   


};

inline void G4ParticleChangeForRadDecay::AddSecondary(G4Track *aTrack)
{
  // add a secondary after size check
  if (theSizeOftheListOfSecondaries > theNumberOfSecondaries) {
    theListOfSecondaries->SetElement(theNumberOfSecondaries, aTrack);
    theNumberOfSecondaries++;
  } else {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << "G4VParticleChange::AddSecondary() Warning  ";
      G4cerr << "theListOfSecondaries is full !! " << G4endl;
      G4cerr << " The object will not be added in theListOfSecondaries" << G4endl;
    }
#endif
  }
}



#endif
















