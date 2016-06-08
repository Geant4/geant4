//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
///////////////////////////////////////////////////////////////////////////////
// File: CCalHit.cc
// Description: Hit class for Calorimeters (Ecal, Hcal, ...)
///////////////////////////////////////////////////////////////////////////////

#include "CCalHit.hh"
#include "g4std/iostream"


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


Hep3Vector   CCalHit::getEntry() const          {return entry;}
void         CCalHit::setEntry(Hep3Vector xyz)  { entry    = xyz; }

double       CCalHit::getIncidentEnergy() const {return theIncidentEnergy; }
void         CCalHit::setIncidentEnergy (double e){theIncidentEnergy  = e; }

int          CCalHit::getTrackID() const         {return theTrackID; }
void         CCalHit::setTrackID (int i)         { theTrackID = i; }

unsigned int CCalHit::getUnitID() const          {return theUnitID; }
void         CCalHit::setUnitID (unsigned int i) { theUnitID = i; }

double       CCalHit::getTimeSlice() const       {return theTimeSlice; }
void         CCalHit::setTimeSlice (double d)    { theTimeSlice = d; }
int          CCalHit::getTimeSliceID() const     {return (int)theTimeSlice;}


void CCalHit::setEnergyDeposit(const double e) { 
  theEnergyDeposit = e; 
}

double CCalHit::getEnergyDeposit() const { 
  return theEnergyDeposit; 
}

void CCalHit::addEnergyDeposit(const CCalHit& aHit) {
  addEnergyDeposit( aHit.getEnergyDeposit() );
}

void  CCalHit::addEnergyDeposit(const double e) { 
  theEnergyDeposit += e; 
}


void CCalHit::print() {
  G4cout << (*this);
}


G4std::ostream& operator<<(G4std::ostream& os, const CCalHit& hit) {
  os << " Data of this CCalHit are:"<< G4endl
     << " \t Time slice ID: " << hit.getTimeSliceID() << G4endl
     << " \t Energy of primary particle (ID = " << hit.getTrackID()
     << ") = " << hit.getIncidentEnergy() << " (MeV)"<< G4endl
     << " \t Entry point in Calorimeter unit number " << hit.getUnitID()
     << " is: " << hit.getEntry() << " (mm)" << G4endl
     << " \t EnergyDeposit = " << hit.getEnergyDeposit() << " (MeV)" << G4endl;
  return os;
}
