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
// $Id: Tst33TimedEventAction.hh,v 1.5 2003-11-18 17:22:30 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// Class Tst33TimedEventAction
//
// Class description:
//
// Tells run manager to stop after time is over. Calculations for
// FOM in special cell.
// 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33TimedEventAction_h
#define Tst33TimedEventAction_h 1

#include "G4Types.hh"
#include "G4Timer.hh"
#include <fstream>

#include "Tst33VEventAction.hh"
#include "Tst33Sigma.hh"

class G4CellScorer;

class Tst33TimedEventAction : public Tst33VEventAction
{
public:
  explicit Tst33TimedEventAction(G4int time);
  virtual ~Tst33TimedEventAction();
  
public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);
  virtual void SpecialCellScorer(const G4CellScorer *scorer);
  virtual void Clear();

private:
  Tst33TimedEventAction(const Tst33TimedEventAction &);
  Tst33TimedEventAction &operator=(const Tst33TimedEventAction &);
  void CalculateFOM();

private:
  G4Timer fTimer;
  const G4CellScorer *fCScorer;
  G4double fProcessTime;
  std::ofstream fOut;
  G4double fMaxRunTime;
  Tst33Sigma fSig;
  G4double fOld_lwe;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
