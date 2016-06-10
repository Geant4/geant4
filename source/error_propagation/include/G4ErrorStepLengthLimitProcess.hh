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
// $Id: G4ErrorStepLengthLimitProcess.hh 66892 2013-01-17 10:57:59Z gunter $
//
//
// Class Description:
//
// Limits the step length if change of direction due to magnetic field
// is too big (user defined limit)

// History:
// - Created:   P. Arce, September 2004
// --------------------------------------------------------------------

#ifndef G4ErrorStepLengthLimitProcess_hh
#define G4ErrorStepLengthLimitProcess_hh

#include "G4ios.hh" 
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ElementTable.hh"
#include "G4Gamma.hh" 
#include "G4Electron.hh"
#include "G4Step.hh" 
#include "G4VErrorLimitProcess.hh"

//-----------------------------------------------------------------

class G4ErrorStepLengthLimitProcess : public G4VErrorLimitProcess
{

 public:  // with description
  
  G4ErrorStepLengthLimitProcess(const G4String& processName =
                                     "G4ErrorStepLengthLimit");
  ~G4ErrorStepLengthLimitProcess();
  
  virtual G4double
  PostStepGetPhysicalInteractionLength( const G4Track& track,
                                              G4double previousStepSize,
                                              G4ForceCondition* condition );
    // returns the step length
};
  
#endif
