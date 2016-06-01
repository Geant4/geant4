// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4FermiConfigurationList_h
#define G4FermiConfigurationList_h 1

#include "globals.hh"
#include "G4FermiConfiguration.hh"
#include "Randomize.hh"

#include <rw/tvvector.h>
#include <rw/tvordvec.h>


class G4FermiConfigurationList 
{
public:
  G4FermiConfigurationList();

  ~G4FermiConfigurationList()
    {};
  
private:
  G4FermiConfigurationList(const G4FermiConfigurationList &right);
  
  const G4FermiConfigurationList & operator=(const G4FermiConfigurationList &right);
  G4bool operator==(const G4FermiConfigurationList &right) const;
  G4bool operator!=(const G4FermiConfigurationList &right) const;
  
public:

  G4bool Initialize(const G4int A, const G4int Z, const G4double TotalEnergyRF);

  G4FermiConfiguration ChooseConfiguration(void);


private:


  enum {MaxNumOfFragments = 6};

  G4double TotNumOfConfigurations; // NumberOfFragments;

  G4double NumOfConfigurations[MaxNumOfFragments]; // NumberOfChannelsPerFragment[MaxNumOfFragments];

  RWTValOrderedVector<G4double> NormalizedWeights;
  
  RWTValOrderedVector<G4FermiConfiguration> Configurations;
};


#endif


