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
// $Id: Tst26EventAction.hh,v 1.3 2003-02-06 11:53:27 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Tst26EventAction_h
#define Tst26EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "g4std/vector"

class Tst26RunAction;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Tst26EventAction : public G4UserEventAction
{
  public:
  
    Tst26EventAction(Tst26RunAction*);
   ~Tst26EventAction();

    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    void SetPrintModulo(G4int    val)  {printModulo = val;};
    void AddEnergy(G4double, G4int, G4int);     

  private:
  
    Tst26RunAction*   Tst26Run;
    G4String          drawFlag;
    G4int             printModulo;          
    G4double          E1;
    G4double          E9;
    G4double          E25;
    G4double          Eabs1;
    G4double          Eabs2;
    G4double          Eabs3;
    G4double          Eabs4;
    G4double          Eth;
    G4std::vector<G4double> Evert;
    G4std::vector<G4int>    Nvert;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
