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
// $Id: Test2Run.hh,v 1.3 2010-09-01 12:56:56 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Test2Run_h
#define Test2Run_h 1

#include "globals.hh"
#include "G4Run.hh"

#include "G4THitsMap.hh"

class G4Event;

class Test2Run : public G4Run
{
public:
  Test2Run();
  virtual ~Test2Run();

public:
  virtual void RecordEvent(const G4Event*);

private:
  G4double GetTotal(const G4THitsMap<G4double> &map) const;

private:
  //-- primitive scorers
  // Maps for accumulation
  // mapSum[7]
  //  i = 0 : eDep
  //  i = 1 : trackLengthGamma
  //  i = 2 : trackLengthElec
  //  i = 3 : trackLengthPosi
  //  i = 4 : nStepGamma
  //  i = 5 : nStepElec
  //  i = 6 : nStepPosi
  G4THitsMap<G4double> mapSum[7];
  G4int colIDSum[7];
  
  //-- traditional sensitive detetector
  G4int fSdID;
  G4double fSdQuantities[7]; // The indices are the same ones as the above primitive scorers

public:
  inline G4double GetTotal(G4int i) const
  { return GetTotal(mapSum[i]); }

  inline G4double GetSdQuantity(G4int i) const {return fSdQuantities[i];}
};

#endif

