// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserLimits.hh,v 1.2 1999-11-16 17:40:50 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
// class G4UserLimits
//
// Class description:
//
// Simple placeholder for user Step limitations

// Author: Paul Kent August 96
// 
// 01-11-97: change GetMaxAllowedStep(), Hisaya Kurashige
// 08-04-98: new data members, mma
//
#ifndef G4USERLIMITS_HH
#define G4USERLIMITS_HH

#include "globals.hh"
class G4Track;

class G4UserLimits
{
public:
  G4UserLimits(G4double ustepMax = DBL_MAX,
               G4double utrakMax = DBL_MAX,
               G4double utimeMax = DBL_MAX,
	       G4double uekinMin = 0.,
	       G4double urangMin = 0. );
  virtual ~G4UserLimits();

public:
  // If a Logical Volume has a G4UserLimits object, 
  //the Step length should be limited as shorter 
  //than MaxAllowedStep in the volume.
  // In the current design, the others limits are irrelavant in tracking
  
  virtual G4double GetMaxAllowedStep(const G4Track&);  
  virtual G4double GetUserMaxTrackLength(const G4Track&) ;
  virtual G4double GetUserMaxTime (const G4Track&);
  virtual G4double GetUserMinEkine(const G4Track&);
  virtual G4double GetUserMinRange(const G4Track&);
  
  virtual void SetMaxAllowedStep(G4double ustepMax);    
  virtual void SetUserMaxTrackLength(G4double utrakMax);
  virtual void SetUserMaxTime(G4double utimeMax);
  virtual void SetUserMinEkine(G4double uekinMin);
  virtual void SetUserMinRange(G4double urangMin);
 
protected:
  G4double fMaxStep;        //max allowed Step size in this volume 
  G4double fMaxTrack;       //max total track length
  G4double fMaxTime;        //max time
  G4double fMinEkine;       //min kinetic energy
  G4double fMinRange;       //min remaining range
};

#include "G4UserLimits.icc"

#endif
