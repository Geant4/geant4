/*
 * G4ITTimeStepper.hh
 *
 *  Created on: 4 juin 2014
 *      Author: kara
 */

#ifndef G4ITTIMESTEPPER_HH_
#define G4ITTIMESTEPPER_HH_

#include "globals.hh"
#include "tls.hh"
#include <map>

class G4VITStepModel;
class G4ITModelHandler;
class G4UserTimeStepAction;
class G4ITGun;
class G4ITTrackingInteractivity;

class G4VScheduler
{
protected:
  G4VScheduler();
  virtual ~G4VScheduler();
private:
  static G4ThreadLocal G4VScheduler* fpInstance;

public:
  static G4VScheduler* Instance();
  virtual void Initialize(){;}
  virtual void Reset(){;}

  virtual void SetVerbose(int){;}

  virtual void SetGun(G4ITGun*){;}

  virtual void Process();

  virtual G4bool IsRunning(){ return false; }

  virtual G4ITModelHandler* GetModelHandler(){ return 0; }

  virtual void RegisterModel(G4VITStepModel*, double){;}

  virtual void SetEndTime(const double){;}

  virtual void SetTimeTolerance(double){;}
  // Two tracks below the time tolerance are supposed to be
  // in the same time slice
  virtual double GetTimeTolerance() const{ return -1;}

  virtual void SetMaxZeroTimeAllowed(int){;}
  virtual int GetMaxZeroTimeAllowed() const{ return -1;}

  virtual void SetTimeSteps(std::map<double, double>*){;}
  virtual void AddTimeStep(double /*startingTime*/, double /*timeStep*/){;}
  virtual void SetDefaultTimeStep(double){;}
  virtual double GetLimitingTimeStep() const {return -1;}
  virtual G4int GetNbSteps() const {return -1;}
  virtual void SetMaxNbSteps(G4int) {;}
  virtual G4int GetMaxNbSteps() const {return 0;}
  virtual G4double GetStartTime() const {return -1;}
  virtual G4double GetEndTime() const {return -1;}
  virtual G4double GetTimeStep() const {return -1;}
  virtual G4double GetPreviousTimeStep() const {return -1;}
  virtual G4double GetGlobalTime() const {return -1;}

  virtual void SetUserAction(G4UserTimeStepAction*) {;}
  virtual G4UserTimeStepAction* GetUserTimeStepAction() const {return 0;}

  virtual void SetInteractivity(G4ITTrackingInteractivity*){;}
  virtual G4ITTrackingInteractivity* GetInteractivity() {return 0;}
};

#endif /* G4ITTIMESTEPPER_HH_ */
