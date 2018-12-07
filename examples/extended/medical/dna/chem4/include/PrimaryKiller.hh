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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//

#ifndef CHEM4_PrimaryKiller_h
#define CHEM4_PrimaryKiller_h 1

#include <G4VPrimitiveScorer.hh>
#include <G4UImessenger.hh>
#include <G4THitsMap.hh>

class G4UIcmdWithADoubleAndUnit;

/** \file PrimaryKiller.hh*/

// Description:
//   Kill the primary particle:
//   - either after a given energy loss
//   - or after the primary particle has reached a given energy

class PrimaryKiller : public G4VPrimitiveScorer,
                      public G4UImessenger
{
private:
  double fELoss; // cumulated energy loss by the primary
  
  double fELossRange_Min; // fELoss from which the primary is killed
  double fELossRange_Max; // fELoss from which the event is aborted
  double fKineticE_Min; // kinetic energy below which the primary is killed

  G4UIcmdWithADoubleAndUnit* fpELossUI;
  G4UIcmdWithADoubleAndUnit* fpAbortEventIfELossUpperThan;
  G4UIcmdWithADoubleAndUnit* fpMinKineticE;
  
public:
  PrimaryKiller(G4String name, G4int depth=0);

  virtual ~PrimaryKiller();
  
  /** Set energy under which the particle should be
   killed*/
  inline void SetEnergyThreshold(double energy){
    fKineticE_Min = energy;
  }
  
  /** Set the energy loss from which the primary is
   killed*/
  inline void SetMinLossEnergyLimit(double energy){
    fELossRange_Min = energy;
  }
  
  /** Set the energy loss from which the event is
   aborted*/
  inline void SetMaxLossEnergyLimit(double energy){
    fELossRange_Max = energy;
  }

  /** Method related to G4UImessenger
      used to control energy cuts through macro file
   */
  virtual void SetNewValue(G4UIcommand * command,
                           G4String newValue);
  
protected:
  virtual G4bool ProcessHits(G4Step*,
                             G4TouchableHistory*);

public:
  virtual void Initialize(G4HCofThisEvent*);
  virtual void EndOfEvent(G4HCofThisEvent*);
  virtual void clear();
  virtual void DrawAll();
  virtual void PrintAll();
};
#endif
