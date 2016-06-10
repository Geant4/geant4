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
// $Id: G4VProcessPlacer.hh 77477 2013-11-25 09:42:24Z gcosmo $
//
// ----------------------------------------------------------------------
// Class G4VProcessPlacer
//
// Class description:
//
// Used internally by importance sampling and scoring to place
// processes as second or last PostStepDoIt processes.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VProcessPlacer_hh
#define G4VProcessPlacer_hh G4VProcessPlacer_hh

#include "G4Types.hh"

class G4VProcess;
 
class G4VProcessPlacer
{

public:  // with description

  G4VProcessPlacer();
  virtual ~G4VProcessPlacer();

  virtual void AddProcessAsLastDoIt(G4VProcess *process) = 0;
    // place a post step do it process such that the 
    // PostStepDoIt function is called last
    // THE ORDER CHANGES BY SUBSEQUENT CALLS     

  virtual void AddProcessAsSecondDoIt(G4VProcess *process) = 0;
    // place a post step do it process such that the 
    // PostStepDoIt function is called second
    // THE ORDER CHANGES BY SUBSEQUENT CALLS         

  virtual void RemoveProcess(G4VProcess *process) = 0;
    // removes a given process 

};

#endif
