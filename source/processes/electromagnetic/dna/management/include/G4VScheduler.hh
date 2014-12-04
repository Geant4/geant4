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
