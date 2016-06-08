// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4FermiBreakUp_h
#define G4FermiBreakUp_h 1

#include "G4VFermiBreakUp.hh"
#include "G4FermiConfiguration.hh"
#include "G4FermiConfigurationList.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

class G4FermiBreakUp : public G4VFermiBreakUp 
{
public:
  G4FermiBreakUp();
  ~G4FermiBreakUp();
  
private:
  G4FermiBreakUp(const G4FermiBreakUp &right);
  
  const G4FermiBreakUp & operator=(const G4FermiBreakUp &right);
  G4bool operator==(const G4FermiBreakUp &right) const;
  G4bool operator!=(const G4FermiBreakUp &right) const;
  
public:
  G4FragmentVector * BreakItUp(const G4Fragment &theNucleus);
};


#endif


