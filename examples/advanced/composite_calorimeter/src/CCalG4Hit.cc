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
// File: CCalG4Hit.cc
// Description: G4 Hit class for Calorimeters (Ecal, Hcal, ...)
///////////////////////////////////////////////////////////////////////////////

#include "CCalG4Hit.hh"
#include <iostream>


CCalG4Hit::CCalG4Hit() : CCalHit(), 
  elem(0.0), hadr(0.0) 
{}


CCalG4Hit::~CCalG4Hit() {}


CCalG4Hit::CCalG4Hit(const CCalG4Hit &right) : 
  G4VHit(right), CCalHit(right),
  elem(right.elem), hadr(right.hadr) 
{}


const CCalG4Hit& CCalG4Hit::operator=(const CCalG4Hit &right) {  
  CCalHit::operator=(right);
  elem  = right.elem;
  hadr  = right.hadr;
  return *this;
}


G4double CCalG4Hit::getEM() const      { return elem; }
void   CCalG4Hit::setEM (G4double e)   { elem = e; }
      
G4double CCalG4Hit::getHadr() const    { return hadr; }
void   CCalG4Hit::setHadr (G4double e) { hadr = e; }
      
void CCalG4Hit::addEnergyDeposit(const CCalG4Hit& aHit) {
  addEnergyDeposit( aHit.getEM(), aHit.getHadr() );
}

void CCalG4Hit::addEnergyDeposit(G4double em, G4double hd) {
  elem  += em; 
  hadr += hd;
  CCalHit::addEnergyDeposit(em+hd);
}


void CCalG4Hit::Print() {
  G4cout << (*this);
}


std::ostream& operator<< (std::ostream& os, const CCalG4Hit& hit) {
  os << static_cast<CCalHit>(hit);
  os << " Data specific of this CCalG4Hit are:" << G4endl
     << " \t EnergyDeposit of EM particles = " << hit.getEM() 
     << " (MeV)" << G4endl
     << " \t EnergyDeposit of HD particles = " << hit.getHadr() 
     << " (MeV)" << G4endl;
  return os;
}

