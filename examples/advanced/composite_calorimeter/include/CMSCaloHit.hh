///////////////////////////////////////////////////////////////////////////////
// File: CMSCaloHit.h
// Date: 06.10.99 Veronique.Lefebure@cern.ch
// 
// Hit class for Calorimeters (Ecal, Hcal, ...)
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
// Modified: 19/02/00 I.G. & D.R. Added output operator (<<)
//                                Added const to get methods
//           20/02/00 I.G. & D.R. addEnergy takes a const CMSCaloHit&
//                                instead of a pointer
///////////////////////////////////////////////////////////////////////////////

#ifndef CMSCaloHit_h
#define CMSCaloHit_h 1

#include <CLHEP/Vector/ThreeVector.h>
#include <iostream>

class CMSCaloHit {
  
public:
  
  CMSCaloHit();
  ~CMSCaloHit();
  CMSCaloHit(const CMSCaloHit &right);
  const CMSCaloHit& operator=(const CMSCaloHit &right);
  int operator==(const CMSCaloHit &right){return 0;}
  
  void print();
  
public:
  
  Hep3Vector   getEntry() const;
  void         setEntry(Hep3Vector xyz);
  
  double       getEM() const;
  void         setEM (double e);
  
  double       getHadr() const;
  void         setHadr (double e);
  
  double       getIncidentEnergy() const;
  void         setIncidentEnergy (double e);
  
  int          getTrackID() const;
  void         setTrackID (int i);
  
  unsigned int getUnitID() const;
  void         setUnitID (unsigned int i);
  
  double       getTimeSlice() const;     
  void         setTimeSlice(double d);
  int          getTimeSliceID() const;     
  
  void         addEnergyDeposit(double em, double hd);
  //void       addEnergyDeposit(CMSCaloHit* aHit);
  void         addEnergyDeposit(const CMSCaloHit& aHit);
  
  double       getEnergyDeposit() const;
  
private:
  
  Hep3Vector   entry;             //Entry point
  double       elem;              //EnergyDeposit of EM particles
  double       hadr;              //EnergyDeposit of HD particles
  double       theIncidentEnergy; //Energy of the primary particle
  int          theTrackID;        //Identification number of the primary
                                 //particle
  unsigned int theUnitID;         //Calorimeter Unit Number
  double       theTimeSlice;      //Time Slice Identification
};

ostream& operator<<(ostream&, const CMSCaloHit&);

#endif
