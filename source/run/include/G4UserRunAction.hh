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
//
// $Id: G4UserRunAction.hh,v 1.7 2003/04/04 16:45:30 asaim Exp $
// GEANT4 tag $Name: geant4-05-01 $
//

#ifndef G4UserRunAction_h
#define G4UserRunAction_h 1

class G4Run;

// class description:
//
//  This is the base class of a user's action class which defines the
// user's action at the begining and the end of each run. The user can
// override the following two methods but the user should not change 
// any of the contents of G4Run object.
//    virtual void BeginOfRunAction(const G4Run* aRun);
//    virtual void EndOfRunAction(const G4Run* aRun);
// The user can override the following method to instanciate his/her own
// concrete Run class. G4Run has a virtual method RecordEvent, so that
// the user can store any information useful to him/her with event statistics.
//    virtual G4Run* GenerateRun();
//  The user's concrete class derived from this class must be set to
// G4RunManager via G4RunManager::SetUserAction() method.
//

class G4UserRunAction
{
  public:
    G4UserRunAction();
    virtual ~G4UserRunAction();

  public:
    virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);
};

#endif


