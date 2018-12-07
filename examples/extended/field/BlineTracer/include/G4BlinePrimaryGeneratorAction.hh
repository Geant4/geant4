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
/// \file field/BlineTracer/include/G4BlinePrimaryGeneratorAction.hh
/// \brief Definition of the G4BlinePrimaryGeneratorAction class
//
//
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
    virtual ~G4BlinePrimaryGeneratorAction();

    virtual  void GeneratePrimaries(G4Event* anEvent);
    inline void SetUserPrimaryAction(G4VUserPrimaryGeneratorAction* anAction)
      { fUserPrimaryAction=anAction; }

  private:

    G4VUserPrimaryGeneratorAction* fUserPrimaryAction;
    G4bool fFirstPartOfBline;
    G4ThreeVector fBlineStartPosition;
    G4double fT0;
};

#endif
