// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Implementation of the HETC88 code into Geant4.
// Evaporation and De-excitation parts
// T. Lampen, Helsinki Institute of Physics, May-2000

#ifndef G4VChargedEvapCh_h
#define G4VChargedEvapCh_h 1

#include "globals.hh"
#include "G4BertiniEvaporationChannel.hh"
#include "Randomize.hh"

class G4BEChargedChannel : public G4BertiniEvaporationChannel
{
public:
  G4BEChargedChannel();
  virtual ~G4BEChargedChannel(); 
  
  virtual void calculateProbability();
  virtual G4DynamicParticle * emit() = 0;
  virtual G4double coulombFactor() = 0;
  G4double coulombFactorForProton();
  G4double qmFactorForProton();
  G4double qmFactorForAlpha();
  G4double sampleKineticEnergy();
  
protected:  
  G4double A;
  G4double spin;
};


#endif
