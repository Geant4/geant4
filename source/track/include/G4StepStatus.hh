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
// $Id: G4StepStatus.hh,v 1.3 2001-07-11 10:08:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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


