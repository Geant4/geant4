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
// $Id: G4StepStatus.hh 68795 2013-04-05 13:24:46Z gcosmo $
//
//
//---------------------------------------------------------------
//
// G4StepStatus.hh
//
// Class Description:
//   This is an enumerator to define possible sources which
//   can define the Step length.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#ifndef G4StepStatus_h
#define G4StepStatus_h 1

//////////////////
enum G4StepStatus
//////////////////
{
  fWorldBoundary,           
    // Step reached the world boundary
  fGeomBoundary,            
    // Step defined by a geometry boundary
  fAtRestDoItProc,          
    // Step defined by a PreStepDoItVector
  fAlongStepDoItProc,       
    // Step defined by a AlongStepDoItVector
  fPostStepDoItProc,        
    // Step defined by a PostStepDoItVector
  fUserDefinedLimit,
    // Step defined by the user Step limit in the logical volume
  fExclusivelyForcedProc,   
    // Step defined by an exclusively forced PostStepDoIt process 
  fUndefined                
    // Step not defined yet
};

#endif


