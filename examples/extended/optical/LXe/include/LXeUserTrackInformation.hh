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
#include "G4VUserTrackInformation.hh"
#include "globals.hh"

#ifndef LXeUserTrackInformation_h
#define LXeUserTrackInformation_h 1

enum LXeTrackStatus { active=1, hitPMT=2, absorbed=4, boundaryAbsorbed=8,
                      hitSphere=16, inactive=14};

/*LXeTrackStatus:
  active: still being tracked
  hitPMT: stopped by being detected in a PMT
  absorbed: stopped by being absorbed with G4OpAbsorbtion
  boundaryAbsorbed: stopped by being aborbed with G4OpAbsorbtion
  hitSphere: track hit the sphere at some point
  inactive: track is stopped for some reason
   -This is the sum of all stopped flags so can be used to remove stopped flags
  
 */

class LXeUserTrackInformation : public G4VUserTrackInformation
{
public:
  LXeUserTrackInformation();
  ~LXeUserTrackInformation();
  
  //Sets the track status to s (does not check validity of flags)
  void SetTrackStatusFlags(int s){status=s;}
  //Does a smart add of track status flags (disabling old flags that conflict)
  //If s conflicts with itself it will not be detected
  void AddTrackStatusFlag(int s);
  
  int GetTrackStatus()const {return status;}
  
  void IncReflections(){reflections++;}
  G4int GetReflectionCount()const {return reflections;}

  void SetForceDrawTrajectory(G4bool b){forcedraw=b;}
  G4bool GetForceDrawTrajectory(){return forcedraw;}

  inline void Print()const{};
private:
  int status;
  G4int reflections;
  G4bool forcedraw;
};

#endif
