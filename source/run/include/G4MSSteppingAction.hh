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
// $Id: G4MSSteppingAction.hh,v 1.1 2006-05-04 19:42:47 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

// class description:
//

//////////////////////
//G4MSSteppingAction
/////////////////////


#ifndef G4MSSteppingAction_h
#define G4MSSteppingAction_h 1

class G4Region;

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4MSSteppingAction : public G4UserSteppingAction
{
  public:
    G4MSSteppingAction();
    virtual ~G4MSSteppingAction();

    void Initialize(G4bool rSens,G4Region* reg);
    virtual void UserSteppingAction(const G4Step*);

  private:
    G4bool regionSensitive;
    G4Region* theRegion;
    G4double length;
    G4double x0;
    G4double lambda;

  public:
    inline G4double GetTotalStepLength() const { return length; }
    inline G4double GetX0() const { return x0; }
    inline G4double GetLambda0() const { return lambda; }
};

#endif
