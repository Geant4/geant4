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
///////////////////////////////////////////////////////////////////////////////
// File: CCalHit.cc
// Description: Hit class for Calorimeters (Ecal, Hcal, ...)
///////////////////////////////////////////////////////////////////////////////

#include "CCalHit.hh"
#include <iostream>


CCalHit::CCalHit() : 
  entry(0), theIncidentEnergy(0.0), theTrackID(-1),
  theUnitID(0), theTimeSlice(0.0), theEnergyDeposit(0.0)
{}


CCalHit::~CCalHit() {}


CCalHit::CCalHit(const CCalHit &right) :
  entry( right.entry ), 
  theIncidentEnergy( right.theIncidentEnergy ),
  theTrackID( right.theTrackID ),
  theUnitID( right.theUnitID ),
  theTimeSlice( right.theTimeSlice ),
  theEnergyDeposit( right.theEnergyDeposit ) 
{}


const CCalHit& CCalHit::operator=(const CCalHit &right) {
  entry = right.entry;
  theIncidentEnergy = right.theIncidentEnergy;
  theTrackID = right.theTrackID;
  theUnitID = right.theUnitID;
  theTimeSlice = right.theTimeSlice;
  theEnergyDeposit = right.theEnergyDeposit;
  return *this;
}


CLHEP::Hep3Vector   CCalHit::getEntry() const          {return entry;}
void         CCalHit::setEntry(CLHEP::Hep3Vector xyz)  { entry    = xyz; }

G4double       CCalHit::getIncidentEnergy() const {return theIncidentEnergy; }
void         CCalHit::setIncidentEnergy (G4double e){theIncidentEnergy  = e; }

G4int          CCalHit::getTrackID() const         {return theTrackID; }
void         CCalHit::setTrackID (G4int i)         { theTrackID = i; }

unsigned int CCalHit::getUnitID() const          {return theUnitID; }
void         CCalHit::setUnitID (unsigned int i) { theUnitID = i; }

G4double       CCalHit::getTimeSlice() const       {return theTimeSlice; }
void         CCalHit::setTimeSlice (G4double d)    { theTimeSlice = d; }
G4int          CCalHit::getTimeSliceID() const     { if ( theTimeSlice > 1.0E9 ) return 999999999;
                                                   return (G4int)theTimeSlice;}

void CCalHit::setEnergyDeposit(const G4double e) { 
  theEnergyDeposit = e; 
}

G4double CCalHit::getEnergyDeposit() const { 
  return theEnergyDeposit; 
}

void CCalHit::addEnergyDeposit(const CCalHit& aHit) {
  addEnergyDeposit( aHit.getEnergyDeposit() );
}

void  CCalHit::addEnergyDeposit(const G4double e) { 
  theEnergyDeposit += e; 
}


void CCalHit::print() {
  G4cout << (*this);
}


std::ostream& operator<<(std::ostream& os, const CCalHit& hit) {
  os << " Data of this CCalHit are:"<< G4endl
     << " \t Time slice ID: " << hit.getTimeSliceID() << G4endl
     << " \t Energy of primary particle (ID = " << hit.getTrackID()
     << ") = " << hit.getIncidentEnergy() << " (MeV)"<< G4endl
     << " \t Entry point in Calorimeter unit number " << hit.getUnitID()
     << " is: " << hit.getEntry() << " (mm)" << G4endl
     << " \t EnergyDeposit = " << hit.getEnergyDeposit() << " (MeV)" << G4endl;
  return os;
}
