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
// $Id: Tst33TimedEventAction.hh,v 1.6 2006-06-29 21:59:57 gunter Exp $
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

    
