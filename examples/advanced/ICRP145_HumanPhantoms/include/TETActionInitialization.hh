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
// TETActionInitialization.hh
// 
// Author: Haegin Han
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Geant4 Contributors: J. Allison and S. Guatelli
//
#ifndef TETActionInitialization_h
#define TETActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "G4String.hh"
#include "TETRunAction.hh"
#include "TETSteppingAction.hh"
#include "TETPrimaryGeneratorAction.hh"

class TETModelImport;

// *********************************************************************
// This UserActionInitialization class initializes UserAction classes
// -- BuildForMaster: instantiate UserRunAction class to be used by
//                    G4MTRunManager. This method is not invoked in
//                    the sequential mode.
// -- Build: instantiate UserPrimaryGeneratorAction, UserRunAction, and
//           UserSteppingAction classes.
// *********************************************************************

class TETActionInitialization : public G4VUserActionInitialization
{
public:
 explicit TETActionInitialization(TETModelImport* tetData,
			           G4String outputFileName);

 virtual ~TETActionInitialization() = default;

 virtual void BuildForMaster() const;
 virtual void Build() const;

private:
 TETModelImport* fTetData;
 G4String fOutput;
};
#endif

    
