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
// $Id: G4ProcessPlacer.hh 77477 2013-11-25 09:42:24Z gcosmo $
//
// ----------------------------------------------------------------------
// Class G4ProcessPlacer
//
// Class description:
//
// Used internally by importance sampling and scoring. 
// See G4VProcessPlacer.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ProcessPlacer_hh
#define G4ProcessPlacer_hh G4ProcessPlacer_hh 

#include "G4Types.hh"
#include "G4String.hh"
#include "G4VProcessPlacer.hh"

class G4ProcessManager;
class G4ProcessVector;

class G4ProcessPlacer : public G4VProcessPlacer
{

public:  // with description

  explicit G4ProcessPlacer(const G4String &particlename);
    // create a process placer for a particle type

  virtual ~G4ProcessPlacer();

  virtual void AddProcessAsLastDoIt(G4VProcess *process);
    // place a post step do it process such that the 
    // PostStepDoIt function is called last
    // THE ORDER CHANGES BY SUBSEQUENT CALLS     

  virtual void AddProcessAsSecondDoIt(G4VProcess *process);
    // place a post step do it process such that the 
    // PostStepDoIt function is called second
    // THE ORDER CHANGES BY SUBSEQUENT CALLS         

  virtual void RemoveProcess(G4VProcess *process);
    // removes a given process 

  enum SecondOrLast
  {
    eSecond = 1,            
    eLast = 0
  };

private:

  G4ProcessManager *GetProcessManager();
  void AddProcessAs(G4VProcess *process, SecondOrLast);

  void PrintProcVec(G4ProcessVector* processVec);
  void PrintAlongStepGPILVec();  
  void PrintAlongStepDoItVec();  
  void PrintPostStepGPILVec();  
  void PrintPostStepDoItVec();  

private:

  G4String fParticleName;

};

#endif
