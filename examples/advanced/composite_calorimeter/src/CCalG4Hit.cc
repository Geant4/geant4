///////////////////////////////////////////////////////////////////////////////
// File: CCalG4Hit.cc
// Description: G4 Hit class for Calorimeters (Ecal, Hcal, ...)
///////////////////////////////////////////////////////////////////////////////

#include "CCalG4Hit.hh"
#include <iostream.h>


CCalG4Hit::CCalG4Hit() : CCalHit(), 
  elem(0.0), hadr(0.0) 
{}


CCalG4Hit::~CCalG4Hit() {}


CCalG4Hit::CCalG4Hit(const CCalG4Hit &right) : CCalHit(right),
  elem(right.elem), hadr(right.hadr) 
{}


const CCalG4Hit& CCalG4Hit::operator=(const CCalG4Hit &right) {  
  CCalHit::operator=(right);
  elem  = right.elem;
  hadr  = right.hadr;
  return *this;
}


double CCalG4Hit::getEM() const      { return elem; }
void   CCalG4Hit::setEM (double e)   { elem = e; }
      
double CCalG4Hit::getHadr() const    { return hadr; }
void   CCalG4Hit::setHadr (double e) { hadr = e; }
      
void CCalG4Hit::addEnergyDeposit(const CCalG4Hit& aHit) {
  addEnergyDeposit( aHit.getEM(), aHit.getHadr() );
}

void CCalG4Hit::addEnergyDeposit(double em, double hd) {
  elem  += em; 
  hadr += hd;
  CCalHit::addEnergyDeposit(em+hd);
}


void CCalG4Hit::Print() {
  cout << (*this);
}


ostream& operator<< (ostream& os, const CCalG4Hit& hit) {
  os << static_cast<CCalHit>(hit);
  os << " Data specific of this CCalG4Hit are:" << endl
     << " \t EnergyDeposit of EM particles = " << hit.getEM() 
     << " (MeV)" << endl
     << " \t EnergyDeposit of HD particles = " << hit.getHadr() 
     << " (MeV)" << endl;
  return os;
}

