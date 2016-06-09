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
// $Id: G4BlinePrimaryGeneratorAction.hh,v 1.1 2003/11/25 09:29:46 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// --------------------------------------------------------------------
//
// Class description:
//
// Defines the primary generator action used for tracing magnetic
// field lines.
// It generates the primary vertex by using the user defined primary
// generator action and replaces only the definition of the particle
// to be tracked by a "charged-geantino".
// The user can define start position for field line tracking in
// his/her own primary generator action.

// --------------------------------------------------------------------
// Author: Laurent Desorgher (desorgher@phim.unibe.ch)
//         Created - 2003-10-06
// --------------------------------------------------------------------
#ifndef G4BlinePrimaryGeneratorAction_h
#define G4BlinePrimaryGeneratorAction_h 1

#include <vector>

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4Event;

class G4BlinePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction 
{
  public:  // with description

    G4BlinePrimaryGeneratorAction();    
    ~G4BlinePrimaryGeneratorAction();

    void GeneratePrimaries(G4Event* anEvent);
    inline void SetUserPrimaryAction(G4VUserPrimaryGeneratorAction* anAction)
      { fUserPrimaryAction=anAction; }

  private:

    G4VUserPrimaryGeneratorAction* fUserPrimaryAction;
    G4bool FirstPartOfBline;
    G4ThreeVector BlineStartPosition;
    G4double T0;
};

#endif
