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
// File: CCalG4Hit.hh
// Description: G4 Hit class for Calorimeters (Ecal, Hcal, ...)
///////////////////////////////////////////////////////////////////////////////

#ifndef CCalG4Hit_h
#define CCalG4Hit_h 1

#include "G4VHit.hh"
#include "CCalHit.hh"


class CCalG4Hit : public G4VHit, public CCalHit {

  friend G4std::ostream& operator<< (G4std::ostream& os, const CCalG4Hit& hit);

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
