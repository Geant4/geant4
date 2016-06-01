// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//



#ifndef G4EvaporationProbability_h
#define G4EvaporationProbability_h 1


#include "G4VEmissionProbability.hh"
#include "G4EvaporationChannel.hh"


class G4EvaporationProbability : public G4VEmissionProbability
{
public:
  // Only available constructor
  G4EvaporationProbability(G4VEvaporationChannel * aChannel) 
    { theChannel = aChannel; };

  ~G4EvaporationProbability() {};
private:  
  // Default constructor
  G4EvaporationProbability() {};

  // Copy constructor
  G4EvaporationProbability(const G4EvaporationProbability &right);

  const G4EvaporationProbability & operator=(const G4EvaporationProbability &right);
  G4bool operator==(const G4EvaporationProbability &right) const;
  G4bool operator!=(const G4EvaporationProbability &right) const;
  
public:
  G4double EmissionProbability(const G4Fragment & fragment, const G4double photonExcitation);

private:

  G4double DostrovskyApproximation(const G4int A, const G4double U);
  G4double BotvinaApproximation(const G4int A, const G4double U);
  G4double NikolaiApproximation(const G4int A, const G4double U);




  G4VEvaporationChannel * theChannel;



};


#endif
