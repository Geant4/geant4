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

#include "g4std/iostream"
#include "CLHEP/Vector/ThreeVector.h"
#include "globals.hh"


class CCalHit {

  friend G4std::ostream& operator<< (G4std::ostream&, const CCalHit&);
  
public:
  
  CCalHit();
  ~CCalHit();
  CCalHit(const CCalHit &right);
  const CCalHit& operator=(const CCalHit &right);
  int operator==(const CCalHit &right){return 0;}
  
  void print();
  
public:
  
  Hep3Vector   getEntry() const;
  void         setEntry(Hep3Vector xyz);
  
  double       getIncidentEnergy() const;
  void         setIncidentEnergy (double e);
  
  int          getTrackID() const;
  void         setTrackID (int i);
  
  unsigned int getUnitID() const;
  void         setUnitID (unsigned int i);
  
  double       getTimeSlice() const;     
  void         setTimeSlice(double d);
  int          getTimeSliceID() const;     
  
  double       getEnergyDeposit() const;
  void         setEnergyDeposit(const double e);
  void         addEnergyDeposit(const CCalHit& aHit);
  void         addEnergyDeposit(const double e);
  
private:
  
  Hep3Vector   entry;             // Entry point
  double       theIncidentEnergy; // Energy of the primary particle
  int          theTrackID;        // Identification number of the primary particle
  unsigned int theUnitID;         // Calorimeter Unit Number
  double       theTimeSlice;      // Time Slice Identification
  double       theEnergyDeposit;  // Cumulated Energy deposit

};

#endif
