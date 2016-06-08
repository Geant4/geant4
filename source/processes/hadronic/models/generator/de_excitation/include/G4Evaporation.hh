// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations

#ifndef G4Evaporation_h
#define G4Evaporation_h 1

#include "globals.hh"
#include "g4rw/tvvector.h"
#include "g4rw/tpordvec.h"


#include "G4ios.hh"
#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4EvaporationChannel.hh"
#include "G4CompetitiveFission.hh"
#include "G4PhotonEvaporation.hh"
#include "G4Fragment.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4NucleiProperties.hh"
#include "Randomize.hh"


class G4Evaporation : public G4VEvaporation
{
public:
  G4Evaporation();
  ~G4Evaporation();

private:
  G4Evaporation(const G4Evaporation &right);

  const G4Evaporation & operator=(const G4Evaporation &right);
  G4bool operator==(const G4Evaporation &right) const;
  G4bool operator!=(const G4Evaporation &right) const;

public:
  G4FragmentVector * BreakItUp(const G4Fragment &theNucleus);

private:

  enum {TotNumberOfChannels = 34,
	NumberOfFissionChannel = TotNumberOfChannels-2,
	NumberOfGammaChannel = TotNumberOfChannels-1,
	NumExcitedStates = 35};



  // Excitation energy levels for each channel
  G4RWTValVector<G4double> ExcitEnergyChann00;  // n
  G4RWTValVector<G4double> ExcitEnergyChann01;  // p
  G4RWTValVector<G4double> ExcitEnergyChann02;  // deuteron
  G4RWTValVector<G4double> ExcitEnergyChann03;  // triton
  G4RWTValVector<G4double> ExcitEnergyChann04;  // He3
  G4RWTValVector<G4double> ExcitEnergyChann05;  // alpha
  G4RWTValVector<G4double> ExcitEnergyChann06;  // He5 
  G4RWTValVector<G4double> ExcitEnergyChann07;  // He6
  G4RWTValVector<G4double> ExcitEnergyChann08;  // Li5
  G4RWTValVector<G4double> ExcitEnergyChann09;  // Li5
  G4RWTValVector<G4double> ExcitEnergyChann10;
  G4RWTValVector<G4double> ExcitEnergyChann11;
  G4RWTValVector<G4double> ExcitEnergyChann12;
  G4RWTValVector<G4double> ExcitEnergyChann13;
  G4RWTValVector<G4double> ExcitEnergyChann14;
  G4RWTValVector<G4double> ExcitEnergyChann15;
  G4RWTValVector<G4double> ExcitEnergyChann16;
  G4RWTValVector<G4double> ExcitEnergyChann17;
  G4RWTValVector<G4double> ExcitEnergyChann18;
  G4RWTValVector<G4double> ExcitEnergyChann19;
  G4RWTValVector<G4double> ExcitEnergyChann20;
  G4RWTValVector<G4double> ExcitEnergyChann21;
  G4RWTValVector<G4double> ExcitEnergyChann22;
  G4RWTValVector<G4double> ExcitEnergyChann23;
  G4RWTValVector<G4double> ExcitEnergyChann24;
  G4RWTValVector<G4double> ExcitEnergyChann25;
  G4RWTValVector<G4double> ExcitEnergyChann26;
  G4RWTValVector<G4double> ExcitEnergyChann27;
  G4RWTValVector<G4double> ExcitEnergyChann28;
  G4RWTValVector<G4double> ExcitEnergyChann29; 
  G4RWTValVector<G4double> ExcitEnergyChann30;
  G4RWTValVector<G4double> ExcitEnergyChann31; 


  // Spin of excitation energy levels for each channel
  G4RWTValVector<G4int> ExcitSpinChann00;
  G4RWTValVector<G4int> ExcitSpinChann01;
  G4RWTValVector<G4int> ExcitSpinChann02;
  G4RWTValVector<G4int> ExcitSpinChann03;
  G4RWTValVector<G4int> ExcitSpinChann04;
  G4RWTValVector<G4int> ExcitSpinChann05;
  G4RWTValVector<G4int> ExcitSpinChann06;
  G4RWTValVector<G4int> ExcitSpinChann07;
  G4RWTValVector<G4int> ExcitSpinChann08;
  G4RWTValVector<G4int> ExcitSpinChann09;
  G4RWTValVector<G4int> ExcitSpinChann10;
  G4RWTValVector<G4int> ExcitSpinChann11;
  G4RWTValVector<G4int> ExcitSpinChann12;
  G4RWTValVector<G4int> ExcitSpinChann13;
  G4RWTValVector<G4int> ExcitSpinChann14;
  G4RWTValVector<G4int> ExcitSpinChann15;
  G4RWTValVector<G4int> ExcitSpinChann16;
  G4RWTValVector<G4int> ExcitSpinChann17;
  G4RWTValVector<G4int> ExcitSpinChann18;
  G4RWTValVector<G4int> ExcitSpinChann19;
  G4RWTValVector<G4int> ExcitSpinChann20;
  G4RWTValVector<G4int> ExcitSpinChann21;
  G4RWTValVector<G4int> ExcitSpinChann22;
  G4RWTValVector<G4int> ExcitSpinChann23;
  G4RWTValVector<G4int> ExcitSpinChann24;
  G4RWTValVector<G4int> ExcitSpinChann25;
  G4RWTValVector<G4int> ExcitSpinChann26;
  G4RWTValVector<G4int> ExcitSpinChann27;
  G4RWTValVector<G4int> ExcitSpinChann28;
  G4RWTValVector<G4int> ExcitSpinChann29; 
  G4RWTValVector<G4int> ExcitSpinChann30;
  G4RWTValVector<G4int> ExcitSpinChann31;


  G4VEvaporationChannel * theChannels[TotNumberOfChannels];


};

#endif





