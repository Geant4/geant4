// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserLimits.hh,v 1.5 2000-02-16 14:08:02 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//
// class G4UserLimits
//
// Class description:
//
// Simple placeholder for user Step limitations
//

// Author: Paul Kent August 96
// 
// 01-11-97: change GetMaxAllowedStep(), Hisaya Kurashige
// 08-04-98: new data members, mma
// 02-16-00: add a member of Type, Hisaya Kurashige 
//
#ifndef G4USERLIMITS_HH
#define G4USERLIMITS_HH

#include "globals.hh"

class G4Track;

class G4UserLimits
{
public:  // with description

  G4UserLimits(G4double ustepMax = DBL_MAX,
               G4double utrakMax = DBL_MAX,
               G4double utimeMax = DBL_MAX,
	       G4double uekinMin = 0.,
	       G4double urangMin = 0. );
  G4UserLimits(const G4String& type = "base",
	       G4double ustepMax = DBL_MAX,
               G4double utrakMax = DBL_MAX,
               G4double utimeMax = DBL_MAX,
	       G4double uekinMin = 0.,
	       G4double urangMin = 0. );
  virtual ~G4UserLimits();

public:  // with description

  // If a Logical Volume has a G4UserLimits object, 
  // the Step length should be limited as shorter 
  // than MaxAllowedStep in the volume.
  // In the current design, the other limits are irrelavant in tracking
  
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

  const G4String & GetType() const;
  void  SetType(const G4String& type);
  // Set/Get type name for UserLimits
  // This type member is supposed to be used to check
  // real class type for each object of UserLimits
  // In other words, users who uses special their UserLimits class
  // derived from this base class should name their class by using
  // this type member.  
 
protected:  // with description

  G4double fMaxStep;        // max allowed Step size in this volume 
  G4double fMaxTrack;       // max total track length
  G4double fMaxTime;        // max time
  G4double fMinEkine;       // min kinetic energy  (only for charged particles)
  G4double fMinRange;       // min remaining range (only for charged particles)

  G4String fType;           // type name
};

#include "G4UserLimits.icc"

#endif
