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
//
// $Id: G4IsotopeProperty.cc 98732 2016-08-09 10:50:57Z gcosmo $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
// **********************************************************************
//      New design using G4VIsotopeTable          5 Oct. 99 H.Kurashige

#include "G4ios.hh"
#include <iomanip>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4IsotopeProperty.hh"
#include "G4DecayTable.hh"

// ######################################################################
// ###                           IsotopeProperty                      ###
// ######################################################################

G4IsotopeProperty::G4IsotopeProperty():
                   fAtomicNumber(0),fAtomicMass(0),
		   fISpin(0),fEnergy(0.0),
		   fLifeTime(-1.0),fDecayTable(0),
		   fMagneticMoment(0.0),
		   fIsomerLevel(-1),
                   fFloatLevelBase(G4Ions::G4FloatLevelBase::no_Float)
{
}


G4IsotopeProperty::~G4IsotopeProperty()
{
  if (fDecayTable != 0) delete fDecayTable;
}

G4IsotopeProperty::G4IsotopeProperty(const  G4IsotopeProperty& right)
{
  fAtomicNumber    = right.fAtomicNumber;
  fAtomicMass      = right.fAtomicMass;
  fISpin           = right.fISpin;
  fMagneticMoment  = right.fMagneticMoment;
  fEnergy          = right.fEnergy;
  fLifeTime        = right.fLifeTime;
  fIsomerLevel     = right.fIsomerLevel;
  fFloatLevelBase  = right.fFloatLevelBase;
  // decay table is not copied because G4DecayTable has no copy constructor
  fDecayTable   = 0;
}

// Assignment operator
G4IsotopeProperty & G4IsotopeProperty::operator=(G4IsotopeProperty& right)
{
  if (this != &right) {
    fAtomicNumber    = right.fAtomicNumber;
    fAtomicMass      = right.fAtomicMass;
    fISpin           = right.fISpin;
    fMagneticMoment  = right.fMagneticMoment;
    fEnergy          = right.fEnergy;
    fLifeTime        = right.fLifeTime;
    fIsomerLevel     = right.fIsomerLevel;
    fFloatLevelBase  = right.fFloatLevelBase;
    // decay table is not copied because G4DecayTable has no copy constructor
    fDecayTable   = 0;
  }
  return *this;
}

 
// equal / unequal operator
G4int G4IsotopeProperty::operator==(const G4IsotopeProperty &right) const
{
  G4bool value = true;
  value = value && ( fAtomicNumber    == right.fAtomicNumber);
  value = value && ( fAtomicMass      == right.fAtomicMass);
  value = value && ( fISpin           == right.fISpin);
  value = value && ( fMagneticMoment  == right.fMagneticMoment);
  value = value && ( fEnergy          == right.fEnergy);
  value = value && ( fLifeTime        == right.fLifeTime);
  value = value && ( fIsomerLevel     == right.fIsomerLevel);
  value = value && ( fFloatLevelBase  == right.fFloatLevelBase);
  return value;
}

G4int G4IsotopeProperty::operator!=(const G4IsotopeProperty &right) const
{
  return !(*this == right);
}

void G4IsotopeProperty::DumpInfo() const
{
#ifdef G4VERBOSE
  G4cout << "AtomicNumber: " << fAtomicNumber << ",  "
	 << "AtomicMass: " << fAtomicMass << G4endl;
  if (fISpin %2){
    G4cout << "Spin: " << fISpin << "/2";
  } else {
    G4cout << "Spin: " << fISpin /2;
  }
  G4cout << ",   " << "MagneticMoment: " 
	 << fMagneticMoment/MeV*tesla << "[MeV/T]" <<G4endl;
  G4cout << "Isomer Level: "
	 << fIsomerLevel
	 << ", Excited Energy: " 
	 << std::setprecision(1) 
	 << fEnergy/keV;
  if(fFloatLevelBase!=G4Ions::G4FloatLevelBase::no_Float)
  { G4cout << " +" << G4Ions::FloatLevelBaseChar(fFloatLevelBase); }
  G4cout << " [keV]" 
	 << ",   "
	 << std::setprecision(6)
	 << "Life Time: " 
	 << fLifeTime/ns << "[ns]"
         << G4endl;
  if (fDecayTable != 0) {
    fDecayTable->DumpInfo();
  } else {
    // G4cout << "Decay Table is not defined !" << G4endl;
  }
#endif
}







