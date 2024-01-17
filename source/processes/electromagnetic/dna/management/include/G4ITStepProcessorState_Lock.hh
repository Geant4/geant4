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
/*
 * G4ITStepProcessorState_Lock.hh
 *
 *  Created on: Nov 6, 2015
 *      Author: mkaramit
 */

#ifndef SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_MANAGEMENT_INCLUDE_G4ITSTEPPROCESSORSTATE_LOCK_HH_
#define SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_MANAGEMENT_INCLUDE_G4ITSTEPPROCESSORSTATE_LOCK_HH_

class G4TrackingInformation;

class G4ITStepProcessorState_Lock
{
  friend class G4TrackingInformation;
protected:
  inline virtual ~G4ITStepProcessorState_Lock()
  = default;
};

#endif /* SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_MANAGEMENT_INCLUDE_G4ITSTEPPROCESSORSTATE_LOCK_HH_ */
