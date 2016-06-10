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
// $Id: G4MSSteppingAction.hh 66892 2013-01-17 10:57:59Z gunter $
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
