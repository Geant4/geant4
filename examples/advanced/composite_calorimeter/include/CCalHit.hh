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
// File: CCalHit.h
// Description: Hit class for Calorimeters (Ecal, Hcal, ...)
//
// One Hit object should be created
// -for each new particle entering the calorimeter
// -for each detector unit (= cristal or fiber or scintillator layer)
// -for each nanosecond of the shower development
//
// This implies that all hit objects created for a given shower
// have the same value for
// - Entry (= local coordinates of the entrance point of the particle
//            in the unit where the shower starts) 
// - the TrackID (= Identification number of the incident particle)
// - the IncidentEnergy (= energy of that particle)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef CCalHit_h
#define CCalHit_h 1

#include <iostream>
#include <CLHEP/Vector/ThreeVector.h>
#include "globals.hh"


class CCalHit {

  friend std::ostream& operator<< (std::ostream&, const CCalHit&);
  
public:
  
  CCalHit();
  ~CCalHit();
  CCalHit(const CCalHit &right);
  const CCalHit& operator=(const CCalHit &right);
  G4bool operator==(const CCalHit &){return false;}
  
  void print();
  
public:
  
  CLHEP::Hep3Vector   getEntry() const;
  void         setEntry(CLHEP::Hep3Vector xyz);
  
  G4double     getIncidentEnergy() const;
  void         setIncidentEnergy (G4double e);
  
  G4int        getTrackID() const;
  void         setTrackID (G4int i);
  
  unsigned int getUnitID() const;
  void         setUnitID (unsigned int i);
  
  G4double       getTimeSlice() const;     
  void         setTimeSlice(G4double d);
  G4int          getTimeSliceID() const;     
  
  G4double     getEnergyDeposit() const;
  void         setEnergyDeposit(const G4double e);
  void         addEnergyDeposit(const CCalHit& aHit);
  void         addEnergyDeposit(const G4double e);
  
private:
  
  CLHEP::Hep3Vector   entry;        // Entry point
  G4double       theIncidentEnergy; // Energy of the primary particle
  G4int          theTrackID;        // Identification number of the primary particle
  unsigned int theUnitID;           // Calorimeter Unit Number
  G4double       theTimeSlice;      // Time Slice Identification
  G4double       theEnergyDeposit;  // Cumulated Energy deposit

};

#endif
