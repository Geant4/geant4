///////////////////////////////////////////////////////////////////////////////
// File: CCalG4Hit.hh
// Description: G4 Hit class for Calorimeters (Ecal, Hcal, ...)
///////////////////////////////////////////////////////////////////////////////

#ifndef CCalG4Hit_h
#define CCalG4Hit_h 1

#include "G4VHit.hh"
#include "CCalHit.hh"


class CCalG4Hit : public G4VHit, public CCalHit {

  friend ostream& operator<< (ostream& os, const CCalG4Hit& hit);

public:

  CCalG4Hit();
  ~CCalG4Hit();
  CCalG4Hit(const CCalG4Hit & right);
  const CCalG4Hit& operator=(const CCalG4Hit & right); 
  int operator==(const CCalG4Hit &right){return 0;}

public:

  void Draw(){}
  void Print();

  double getEM() const;
  void   setEM (double e);
  
  double getHadr() const;
  void   setHadr (double e);
  
  void  addEnergyDeposit(const CCalG4Hit& aHit);
  void  addEnergyDeposit(double em, double hd);

private:

  double elem; // EnergyDeposit of EM particles
  double hadr; // EnergyDeposit of HD particles

};
#endif
