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
// $Id: ExN07Run.hh,v 1.5 2006/06/29 17:54:50 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//

#ifndef ExN07Run_h
#define ExN07Run_h 1

#include "globals.hh"
#include "G4Run.hh"

#include "G4THitsMap.hh"

class G4Event;

class ExN07Run : public G4Run
{
  public:
    ExN07Run();
    virtual ~ExN07Run();

  public:
    virtual void RecordEvent(const G4Event*);

  private:
    G4double GetTotal(const G4THitsMap<G4double> &map) const;
    G4double FindMinimum(const G4THitsMap<G4double> &map) const;

  private:
    // Maps for accumulation
    // mapSum[i][j]
    //  i = 0 : Calor-A_abs    j = 0 : total eDep
    //  i = 1 : Calor-A_gap    j = 1 : number of gamma
    //  i = 2 : Calor-B_abs    j = 2 : number of electron
    //  i = 3 : Calor-B_gap    j = 3 : number of positron
    //  i = 4 : Calor-C_abs    j = 4 : total step length for e+/e-
    //  i = 5 : Calor-C_gap    j = 5 : total number of steps for e+/e-
    G4THitsMap<G4double> mapSum[6][6];
    G4int colIDSum[6][6];

    // Maps for minimum value
    //  i = 0 : Calor-A_abs    j = 0 : minimum kinE at generation for gamma
    //  i = 1 : Calor-A_gap    j = 1 : minimum kinE at generation for electron
    //  i = 2 : Calor-B_abs    j = 2 : minimum kinE at generation for positron
    //  i = 3 : Calor-B_gap
    //  i = 4 : Calor-C_abs
    //  i = 5 : Calor-C_gap
    G4THitsMap<G4double> mapMin[6][3];
    G4int colIDMin[6][3];

  public:
    inline G4double GetTotalE(G4int i) const
    { return GetTotal(mapSum[i][0]); }
    inline G4double GetNGamma(G4int i) const
    { return GetTotal(mapSum[i][1]); }
    inline G4double GetNElectron(G4int i) const
    { return GetTotal(mapSum[i][2]); }
    inline G4double GetNPositron(G4int i) const
    { return GetTotal(mapSum[i][3]); }
    inline G4double GetTotalL(G4int i) const
    { return GetTotal(mapSum[i][4]); }
    inline G4double GetNStep(G4int i) const
    { return GetTotal(mapSum[i][5]); }
    inline G4double GetEMinGamma(G4int i) const
    { return FindMinimum(mapMin[i][0]); }
    inline G4double GetEMinElectron(G4int i) const
    { return FindMinimum(mapMin[i][1]); }
    inline G4double GetEMinPositron(G4int i) const
    { return FindMinimum(mapMin[i][2]); }
};

#endif

