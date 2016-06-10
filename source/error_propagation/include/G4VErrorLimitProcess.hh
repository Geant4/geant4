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
// $Id: G4VErrorLimitProcess.hh 66892 2013-01-17 10:57:59Z gunter $
//
//
// Class description:
//
// Base class for classes limiting the step.

// History:
//
// Created:     P.Arce          May 2007
// --------------------------------------------------------------------

#ifndef G4VErrorLimitProcess_hh
#define G4VErrorLimitProcess_hh

#include "G4ios.hh" 
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ElementTable.hh"
#include "G4Gamma.hh" 
#include "G4Electron.hh"
#include "G4Step.hh" 

class G4ErrorLimitsMessenger;

//-----------------------------------------------------------------
 
class G4VErrorLimitProcess : public G4VDiscreteProcess
{

 public:  // with description
  
  G4VErrorLimitProcess(const G4String& processName);
  
  ~G4VErrorLimitProcess();
  
  virtual G4double
  PostStepGetPhysicalInteractionLength( const G4Track& track,
                                              G4double previousStepSize,
                                              G4ForceCondition* condition ) = 0;
    // Returns the step limit

  virtual  G4double GetMeanFreePath(const class G4Track &, G4double,
                                    enum G4ForceCondition *);
    // Returns kInfinity

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& );
    // No action but retrieving the G4VParticleChange
    // extracted from the G4Track

  // Get and Set methods

  G4double GetStepLimit() const { return theStepLimit; }
  void SetStepLimit( G4double val ) { theStepLimit = val; }

 protected:

  G4double theStepLimit;  // limit set by the user
  G4double theStepLength; // step length extracted from the user step limit,
                          // with the algorithms of each concrete class

  G4VParticleChange theParticleChange;
};
  
#endif
