// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPDataPoint.hh,v 1.2 1999-06-29 18:43:51 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPDataPoint_h
#define G4NeutronHPDataPoint_h 1

#include "globals.hh"

class G4NeutronHPDataPoint
{
  public:
  
  G4NeutronHPDataPoint();
  G4NeutronHPDataPoint(G4double e, G4double x);
  
  void operator= (const G4NeutronHPDataPoint & aSet)
  {
    if(&aSet!=this)
    {
      energy = aSet.GetEnergy();
      xSec   = aSet.GetXsection();
    }
  }

//  ~G4NeutronHPDataPoint(){}
  
  inline G4double GetEnergy() const   {return energy;}
  inline G4double GetXsection() const {return xSec;}
  
  inline void SetEnergy(G4double e)  {energy = e;}
  inline void SetXsection(G4double x){xSec = x;}
  
  inline G4double GetX() const {return energy;}
  inline G4double GetY() const {return xSec;}
  
  inline void SetX(G4double e)  {energy = e;}
  inline void SetY(G4double x)  {xSec = x;}
  
  inline void SetData(G4double e, G4double x) {energy = e; xSec = x;}
  
  private:
  
  G4double energy;
  G4double xSec;
};

#endif
