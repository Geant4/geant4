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
