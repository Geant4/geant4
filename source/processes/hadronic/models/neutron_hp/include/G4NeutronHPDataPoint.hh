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
//
// $Id: G4NeutronHPDataPoint.hh,v 1.4 2001-07-11 10:06:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPDataPoint_h
#define G4NeutronHPDataPoint_h 1

#include "globals.hh"

class G4NeutronHPDataPoint
{
  public:
  
  G4NeutronHPDataPoint(){energy = 0; xSec = 0;}
  G4NeutronHPDataPoint(G4double e, G4double x){ energy = e; xSec = x;}
  
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
