// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IsotopeProperty.cc,v 1.1 1999-10-05 06:45:19 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
// **********************************************************************
//      New design using G4VIsotopeTable          5 Oct. 99 H.Kurashige

#include "G4ios.hh"
#include <iomanip.h>

#include "G4IsotopeProperty.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                           IsotopeProperty                      ###
// ######################################################################

G4IsotopeProperty::G4IsotopeProperty():
                   fAtomicNumber(0),fAtomicMass(0),
		   fISpin(0),fEnergy(0.0),
		   fLifeTime(-1.0),fDecayTable(0)
{
}


G4IsotopeProperty::~G4IsotopeProperty()
{
  if (fDecayTable != 0) delete fDecayTable;
}

G4IsotopeProperty::G4IsotopeProperty(const  G4IsotopeProperty& right)
{
  fAtomicNumber = right.fAtomicNumber;
  fAtomicMass   = right.fAtomicMass;
  fISpin        = right.fISpin;
  fEnergy       = right.fEnergy;
  fLifeTime     = right.fLifeTime;
  // decay table is not copied because G4DecayTable has no copy constructor
  fDecayTable   = 0;
}

// Assignment operator
G4IsotopeProperty & G4IsotopeProperty::operator=(G4IsotopeProperty& right)
{
  if (this != &right) {
    fAtomicNumber = right.fAtomicNumber;
    fAtomicMass   = right.fAtomicMass;
    fISpin        = right.fISpin;
    fEnergy       = right.fEnergy;
    fLifeTime     = right.fLifeTime;
    // decay table is not copied because G4DecayTable has no copy constructor
    fDecayTable   = 0;
  }
  return *this;
}

 
// equal / unequal operator
G4int G4IsotopeProperty::operator==(const G4IsotopeProperty &right) const
{
  G4bool value = true;
  value = value && ( fAtomicNumber == right.fAtomicNumber);
  value = value && ( fAtomicMass   == right.fAtomicMass);
  value = value && ( fISpin        == right.fISpin);
  value = value && ( fEnergy       == right.fEnergy);
  value = value && ( fLifeTime     == right.fLifeTime);
  return value;
}
G4int G4IsotopeProperty::operator!=(const G4IsotopeProperty &right) const
{
  return !(*this == right);
}

void G4IsotopeProperty::DumpInfo() const
{
  G4cout << "AtomicNumber: " << fAtomicNumber << endl;
  G4cout << "AtomicMass: " << fAtomicMass << endl;
  G4cout << "Spin: " << fISpin << "/2" << endl;
  G4cout << "Excited Energy: " << setprecision(1) << fEnergy/keV << "[keV]" << endl;
  G4cout << "Life Time: " << fLifeTime/ns << "[ns]" << endl;
  if (fDecayTable != 0) {
    fDecayTable->DumpInfo();
  } else {
    G4cout << "Decay Table is not defined !" << endl;
  }
}







