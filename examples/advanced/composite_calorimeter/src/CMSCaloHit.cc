///////////////////////////////////////////////////////////////////////////////
// File: CMSCaloHit.cc
// Date: 18/06/99 Mikhail Kossov
// Modified: 18/06/99 M.Kossov
//           13/07/99 S.B.  Change of name
//           19/02/00 I.G. & D.R. Added output operator (<<)
//                                Added const to get methods
///////////////////////////////////////////////////////////////////////////////

#include "CMSCaloHit.hh"
#include <iostream.h>

CMSCaloHit::CMSCaloHit():entry(0) {

  elem     = 0.;
  hadr     = 0.;
  theIncidentEnergy = 0.;
  theTrackID = -1;
  theUnitID  =  0;
  theTimeSlice = 0.;
}


CMSCaloHit::~CMSCaloHit(){}


CMSCaloHit::CMSCaloHit(const CMSCaloHit &right) {

  elem     = right.elem;
  hadr     = right.hadr;
  theIncidentEnergy  = right.theIncidentEnergy;
  theTrackID = right.theTrackID;
  theUnitID = right.theUnitID;
  theTimeSlice = right.theTimeSlice;
  entry    = right.entry;
}


const CMSCaloHit& CMSCaloHit::operator=(const CMSCaloHit &right) {

  elem     = right.elem;
  hadr     = right.hadr;
  theIncidentEnergy = right.theIncidentEnergy;
  theTrackID = right.theTrackID;
  theUnitID = right.theUnitID;
  theTimeSlice = right.theTimeSlice;
  entry    = right.entry;
 
  return *this;
}


void CMSCaloHit::addEnergyDeposit(const CMSCaloHit& aHit) {

  elem += aHit.getEM();
  hadr += aHit.getHadr();
}


void CMSCaloHit::print() {
  cout << (*this);
}


Hep3Vector   CMSCaloHit::getEntry() const          {return entry;}
void         CMSCaloHit::setEntry(Hep3Vector xyz)  { entry    = xyz; }

double       CMSCaloHit::getEM() const             {return elem; }
void         CMSCaloHit::setEM (double e)          { elem     = e; }
      
double       CMSCaloHit::getHadr() const           {return hadr; }
void         CMSCaloHit::setHadr (double e)        { hadr     = e; }
      
double       CMSCaloHit::getIncidentEnergy() const {return theIncidentEnergy; }
void         CMSCaloHit::setIncidentEnergy (double e){theIncidentEnergy  = e; }

int          CMSCaloHit::getTrackID() const         {return theTrackID; }
void         CMSCaloHit::setTrackID (int i)         { theTrackID = i; }

unsigned int CMSCaloHit::getUnitID() const          {return theUnitID; }
void         CMSCaloHit::setUnitID (unsigned int i) { theUnitID = i; }

double       CMSCaloHit::getTimeSlice() const       {return theTimeSlice; }
void         CMSCaloHit::setTimeSlice (double d)    { theTimeSlice = d; }
int          CMSCaloHit::getTimeSliceID() const     {return (int)theTimeSlice;}

void         CMSCaloHit::addEnergyDeposit(double em, double hd)
                                           {elem  += em ; hadr += hd;}

double       CMSCaloHit::getEnergyDeposit() const      {return elem+hadr;}

ostream& operator<<(ostream& os, const CMSCaloHit& hit) {
  os << " Data of this CMSCaloHit are:"<<endl
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


