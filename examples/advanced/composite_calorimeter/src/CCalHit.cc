///////////////////////////////////////////////////////////////////////////////
// File: CCalHit.cc
// Description: Hit class for Calorimeters (Ecal, Hcal, ...)
///////////////////////////////////////////////////////////////////////////////

#include "CCalHit.hh"
#include <iostream.h>

CCalHit::CCalHit():entry(0) {

  elem     = 0.;
  hadr     = 0.;
  theIncidentEnergy = 0.;
  theTrackID = -1;
  theUnitID  =  0;
  theTimeSlice = 0.;
}


CCalHit::~CCalHit(){}


CCalHit::CCalHit(const CCalHit &right) {

  elem     = right.elem;
  hadr     = right.hadr;
  theIncidentEnergy  = right.theIncidentEnergy;
  theTrackID = right.theTrackID;
  theUnitID = right.theUnitID;
  theTimeSlice = right.theTimeSlice;
  entry    = right.entry;
}


const CCalHit& CCalHit::operator=(const CCalHit &right) {

  elem     = right.elem;
  hadr     = right.hadr;
  theIncidentEnergy = right.theIncidentEnergy;
  theTrackID = right.theTrackID;
  theUnitID = right.theUnitID;
  theTimeSlice = right.theTimeSlice;
  entry    = right.entry;
 
  return *this;
}


void CCalHit::addEnergyDeposit(const CCalHit& aHit) {

  elem += aHit.getEM();
  hadr += aHit.getHadr();
}


void CCalHit::print() {
  cout << (*this);
}


Hep3Vector   CCalHit::getEntry() const          {return entry;}
void         CCalHit::setEntry(Hep3Vector xyz)  { entry    = xyz; }

double       CCalHit::getEM() const             {return elem; }
void         CCalHit::setEM (double e)          { elem     = e; }
      
double       CCalHit::getHadr() const           {return hadr; }
void         CCalHit::setHadr (double e)        { hadr     = e; }
      
double       CCalHit::getIncidentEnergy() const {return theIncidentEnergy; }
void         CCalHit::setIncidentEnergy (double e){theIncidentEnergy  = e; }

int          CCalHit::getTrackID() const         {return theTrackID; }
void         CCalHit::setTrackID (int i)         { theTrackID = i; }

unsigned int CCalHit::getUnitID() const          {return theUnitID; }
void         CCalHit::setUnitID (unsigned int i) { theUnitID = i; }

double       CCalHit::getTimeSlice() const       {return theTimeSlice; }
void         CCalHit::setTimeSlice (double d)    { theTimeSlice = d; }
int          CCalHit::getTimeSliceID() const     {return (int)theTimeSlice;}

void         CCalHit::addEnergyDeposit(double em, double hd)
                                           {elem  += em ; hadr += hd;}

double       CCalHit::getEnergyDeposit() const      {return elem+hadr;}

ostream& operator<<(ostream& os, const CCalHit& hit) {
  os << " Data of this CCalHit are:"<<endl
     << " Time slice ID: " << hit.getTimeSliceID() << endl
     << " EnergyDeposit of EM particles = " << hit.getEM() << endl
     << " EnergyDeposit of HD particles = " << hit.getHadr() << endl
     << " Energy of primary particle (ID = " << hit.getTrackID()
     << ") = " << hit.getIncidentEnergy() << " (MeV)"<<endl
     << " Entry point in Calorimeter unit number " << hit.getUnitID()
     << " is: " << hit.getEntry() << " (mm)" << endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
     << endl;
  return os;
}


