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
// $Id: G4UserRunAction.hh 95232 2016-02-01 14:31:22Z gcosmo $
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
#include "G4Types.hh"

class G4UserRunAction
{
  public:
    G4UserRunAction();
    virtual ~G4UserRunAction();

  public:
    virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

  protected:
    G4bool isMaster;

  public:
    inline virtual void SetMaster(G4bool val=true)
    { isMaster = val; }
    inline G4bool IsMaster() const
    { return isMaster; }
};

#endif


