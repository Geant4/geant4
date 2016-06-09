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
// $Id: RemSimRunAction.hh,v 1.7 2004/05/27 10:13:54 guatelli Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//
//    ********************************************
//    *                                          *
//    *      RemSimRunAction.hh                  *
//    *                                          *
//    ********************************************
//
#ifndef RemSimRunAction_h
#define RemSimRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4RunManager.hh"
#include "globals.hh"

class G4Run;
class RemSimAnalysisManager;
class RemSimRunMessenger;

class RemSimRunAction : public G4UserRunAction
{
public:
  RemSimRunAction();
  ~RemSimRunAction();

public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run* );
  void Read(G4String);
  void ReadData(G4double, G4String);

  G4DataVector* GetPrimaryParticleEnergy();
  G4DataVector* GetPrimaryParticleEnergyDistribution();
  G4double GetPrimaryParticleEnergyDistributionSum();
  G4int GetRunID() {return runID;};
  G4bool GetFile();

private:
  G4DataVector* energies;
  G4DataVector* data;
  G4int runID;
  RemSimRunMessenger* messenger;
  G4String file;
};
#endif



