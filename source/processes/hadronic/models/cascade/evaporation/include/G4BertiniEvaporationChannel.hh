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

#ifndef G4VEvaporationChannel_h
#define G4VEvaporationChannel_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4NucleiProperties.hh" 

class G4BertiniEvaporationChannel
{
public:
  G4BertiniEvaporationChannel();
  virtual ~G4BertiniEvaporationChannel(); 

  void setVerboseLevel( G4int verbose );
  void setNucleusA( G4int inputA );
  void setNucleusZ( G4int inputZ );
  G4int getParticleA();
  G4int getParticleZ();
  void setExcitationEnergy( G4double inputE );
  void setQ( G4double inputQ );
  void setPairingCorrection( G4int isCorrection );
  G4double getLevelDensityParameter();
  G4String getName();
  
  virtual G4double getProbability() ;
  virtual void setProbability( G4double newProb ) ;
  virtual void calculateProbability() = 0 ;
  virtual G4double qmFactor ();
  virtual G4double getQ();
  virtual G4double getCoulomb();
  virtual G4double getThresh();

  virtual G4DynamicParticle* emit() = 0;
  
protected:  
  G4String name;
  G4int verboseLevel;
  G4int nucleusA;
  G4int nucleusZ;
  G4int particleA;
  G4int particleZ;
  G4double exmass;
  G4double emissionProbability;
  G4double rho;
  G4double correction;
  G4double excitationEnergy;
  //  G4int massInNeutronMasses;
  G4int spin;
  G4double Q( G4double a, G4double z );
  G4double pairingEnergyProtons( G4int A );
  G4double pairingEnergyNeutrons( G4int N );
  G4double cameron( G4double a, G4double z );
  G4double cameronShellCorrectionP( G4int Z );
  G4double cameronShellCorrectionN( G4int N );
  void isotropicCosines( G4double&,
			 G4double&,
			 G4double& );
};


#endif
