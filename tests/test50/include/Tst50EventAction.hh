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
// $Id: Tst50EventAction.hh,v 1.10 2003-04-28 14:58:56 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 
#ifndef Tst50EventAction_h
#define Tst50EventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"
class Tst50TrackerHit;
class Tst50AnalysisManager;
class G4Event;
class Tst50PrimaryGeneratorAction;


class Tst50EventAction : public G4UserEventAction
{
public:
  Tst50EventAction(Tst50PrimaryGeneratorAction*);
  ~Tst50EventAction();

  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  G4int    GetEventno();
  G4double RadiationYield(G4double);
  void CalculateEnergyDeposit(G4double);
  void CalculateEnergyDepositPrimary(G4double);
  void CalculateEnergyDepositSecondary(G4double);
  void Number_Steps();
  void Track_length(G4double);

private:
  G4int hit_CollID;
  G4double energy;
  G4double energyDep;
  Tst50PrimaryGeneratorAction* p_Primary;
  G4String filename;
  G4bool RadiationY; 
  G4double energyDepPrimary;
  G4double energyDepSecondary;
  G4bool Foil;
  G4int N_Steps;
  G4double length;
};

#endif

    
