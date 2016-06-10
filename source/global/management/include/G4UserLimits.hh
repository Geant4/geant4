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
//
// $Id: G4UserLimits.hh 67970 2013-03-13 10:10:06Z gcosmo $
// 
//
// class G4UserLimits
//
// Class description:
//
// Simple placeholder for user Step limitations
// In order to activate these limitations, users need to register
// their "special" processes to each particle they want. 
// Sample processes below can be found under processes/transportation 
//   MaxAllowedStep    : UserStepLimit
//   other limitation  : UserSpecialCuts
// In addition, users can add their own Step limitations by creating 
// new class derived from G4UserLimits. In these case, fType member
// is supposed to be used to identify class. 
//  
// Author: Paul Kent August 96
// 
// 01-11-97: change GetMaxAllowedStep(), Hisaya Kurashige
// 08-04-98: new data members, mma
// 02-16-00: add fType member and accessors, Hisaya Kurashige 
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
  G4UserLimits(const G4String& type,
	       G4double ustepMax = DBL_MAX,
               G4double utrakMax = DBL_MAX,
               G4double utimeMax = DBL_MAX,
	       G4double uekinMin = 0.,
	       G4double urangMin = 0. );
  virtual ~G4UserLimits();

public:  // with description
  
  virtual G4double GetMaxAllowedStep(const G4Track&);  
  // If a Logical Volume has a G4UserLimits object, the Step length can
  // be limited as shorter than MaxAllowedStep in the volume.
 
  // 
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
  // Set/Get type name for UserLimits.
  // This type member is supposed to be used to check real class types for
  // each concrete instantiation of G4UserLimits. In other words, users who
  // use special classes derived from this base class should name their class
  // with a proper identifier.  
 
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
