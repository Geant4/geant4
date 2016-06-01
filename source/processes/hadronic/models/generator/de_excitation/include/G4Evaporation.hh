// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations

#ifndef G4Evaporation_h
#define G4Evaporation_h 1

#include "globals.hh"
#include <rw/tvvector.h>
#include <rw/tpordvec.h>


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
  RWTValVector<G4double> ExcitEnergyChann00;  // n
  RWTValVector<G4double> ExcitEnergyChann01;  // p
  RWTValVector<G4double> ExcitEnergyChann02;  // deuteron
  RWTValVector<G4double> ExcitEnergyChann03;  // triton
  RWTValVector<G4double> ExcitEnergyChann04;  // He3
  RWTValVector<G4double> ExcitEnergyChann05;  // alpha
  RWTValVector<G4double> ExcitEnergyChann06;  // He5 
  RWTValVector<G4double> ExcitEnergyChann07;  // He6
  RWTValVector<G4double> ExcitEnergyChann08;  // Li5
  RWTValVector<G4double> ExcitEnergyChann09;  // Li5
  RWTValVector<G4double> ExcitEnergyChann10;
  RWTValVector<G4double> ExcitEnergyChann11;
  RWTValVector<G4double> ExcitEnergyChann12;
  RWTValVector<G4double> ExcitEnergyChann13;
  RWTValVector<G4double> ExcitEnergyChann14;
  RWTValVector<G4double> ExcitEnergyChann15;
  RWTValVector<G4double> ExcitEnergyChann16;
  RWTValVector<G4double> ExcitEnergyChann17;
  RWTValVector<G4double> ExcitEnergyChann18;
  RWTValVector<G4double> ExcitEnergyChann19;
  RWTValVector<G4double> ExcitEnergyChann20;
  RWTValVector<G4double> ExcitEnergyChann21;
  RWTValVector<G4double> ExcitEnergyChann22;
  RWTValVector<G4double> ExcitEnergyChann23;
  RWTValVector<G4double> ExcitEnergyChann24;
  RWTValVector<G4double> ExcitEnergyChann25;
  RWTValVector<G4double> ExcitEnergyChann26;
  RWTValVector<G4double> ExcitEnergyChann27;
  RWTValVector<G4double> ExcitEnergyChann28;
  RWTValVector<G4double> ExcitEnergyChann29; 
  RWTValVector<G4double> ExcitEnergyChann30;
  RWTValVector<G4double> ExcitEnergyChann31; 


  // Spin of excitation energy levels for each channel
  RWTValVector<G4int> ExcitSpinChann00;
  RWTValVector<G4int> ExcitSpinChann01;
  RWTValVector<G4int> ExcitSpinChann02;
  RWTValVector<G4int> ExcitSpinChann03;
  RWTValVector<G4int> ExcitSpinChann04;
  RWTValVector<G4int> ExcitSpinChann05;
  RWTValVector<G4int> ExcitSpinChann06;
  RWTValVector<G4int> ExcitSpinChann07;
  RWTValVector<G4int> ExcitSpinChann08;
  RWTValVector<G4int> ExcitSpinChann09;
  RWTValVector<G4int> ExcitSpinChann10;
  RWTValVector<G4int> ExcitSpinChann11;
  RWTValVector<G4int> ExcitSpinChann12;
  RWTValVector<G4int> ExcitSpinChann13;
  RWTValVector<G4int> ExcitSpinChann14;
  RWTValVector<G4int> ExcitSpinChann15;
  RWTValVector<G4int> ExcitSpinChann16;
  RWTValVector<G4int> ExcitSpinChann17;
  RWTValVector<G4int> ExcitSpinChann18;
  RWTValVector<G4int> ExcitSpinChann19;
  RWTValVector<G4int> ExcitSpinChann20;
  RWTValVector<G4int> ExcitSpinChann21;
  RWTValVector<G4int> ExcitSpinChann22;
  RWTValVector<G4int> ExcitSpinChann23;
  RWTValVector<G4int> ExcitSpinChann24;
  RWTValVector<G4int> ExcitSpinChann25;
  RWTValVector<G4int> ExcitSpinChann26;
  RWTValVector<G4int> ExcitSpinChann27;
  RWTValVector<G4int> ExcitSpinChann28;
  RWTValVector<G4int> ExcitSpinChann29; 
  RWTValVector<G4int> ExcitSpinChann30;
  RWTValVector<G4int> ExcitSpinChann31;


  G4VEvaporationChannel * theChannels[TotNumberOfChannels];


};

#endif





