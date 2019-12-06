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
// File: CCalG4Hit.hh
// Description: G4 Hit class for Calorimeters (Ecal, Hcal, ...)
///////////////////////////////////////////////////////////////////////////////

#ifndef CCalG4Hit_h
#define CCalG4Hit_h 1

#include "G4VHit.hh"
#include "CCalHit.hh"


class CCalG4Hit : public G4VHit, public CCalHit
{

  friend std::ostream& operator<< (std::ostream& os, const CCalG4Hit& hit);

public:

  CCalG4Hit();
  ~CCalG4Hit();
  CCalG4Hit(const CCalG4Hit & right);
  const CCalG4Hit& operator=(const CCalG4Hit & right); 
  G4bool operator==(const CCalG4Hit &){return 0;}

public:

  void Draw(){}
  void Print();

  G4double getEM() const;
  void   setEM (double e);
  
  G4double getHadr() const;
  void   setHadr (G4double e);
  
  void  addEnergyDeposit(const CCalG4Hit& aHit);
  void  addEnergyDeposit(G4double em, G4double hd);

private:

  G4double elem; // EnergyDeposit of EM particles
  G4double hadr; // EnergyDeposit of HD particles

};
#endif
