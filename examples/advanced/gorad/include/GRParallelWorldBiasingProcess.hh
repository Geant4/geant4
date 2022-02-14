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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRParallelWorldBiasingProcess.hh
//   Header file of a process that takes care of the geometry
//   importance biasing
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#ifndef GRParallelWorldBiasingProcess_h
#define GRParallelWorldBiasingProcess_h 1

#include "globals.hh"
class G4Step;
class G4Track;
class G4VParticleChange;
class G4ParticleChange;
class G4ParticleChangeForNothing;
#include "G4ParallelWorldProcess.hh"

// Class Description:

class GRParallelWorldBiasingProcess : public G4ParallelWorldProcess
{
public: 
  GRParallelWorldBiasingProcess(const G4String& processName = "ParaWorld",
                                     G4ProcessType theType = fParallel);
  virtual ~GRParallelWorldBiasingProcess();
  
  G4VParticleChange* PostStepDoIt(const G4Track&,const G4Step&);

private:
  G4ParticleChange* particleChange;
  G4ParticleChangeForNothing* emptyParticleChange;

};

#endif
