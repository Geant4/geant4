///////////////////////////////////////////////////////////////////////////////
// File: CCalHit.cc
// Description: Hit class for Calorimeters (Ecal, Hcal, ...)
///////////////////////////////////////////////////////////////////////////////

#include "CCalHit.hh"
#include <iostream.h>


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
  cout << (*this);
}


ostream& operator<<(ostream& os, const CCalHit& hit) {
  os << " Data of this CCalHit are:"<<endl
     << " \t Time slice ID: " << hit.getTimeSliceID() << endl
     << " \t Energy of primary particle (ID = " << hit.getTrackID()
     << ") = " << hit.getIncidentEnergy() << " (MeV)"<<endl
     << " \t Entry point in Calorimeter unit number " << hit.getUnitID()
     << " is: " << hit.getEntry() << " (mm)" << endl
     << " \t EnergyDeposit = " << hit.getEnergyDeposit() << " (MeV)" << endl;
  return os;
}
