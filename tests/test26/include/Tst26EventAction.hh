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
// $Id: Tst26EventAction.hh,v 1.5 2003-06-16 17:15:45 gunter Exp $
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
#include <vector>

class Tst26RunAction;
class Tst26EventMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Tst26EventAction : public G4UserEventAction
{
  public:
  
    Tst26EventAction(Tst26RunAction*);
   ~Tst26EventAction();

    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    void SetPrintModulo(G4int    val)  {printModulo = val;};
    void SetDrawFlag   (G4String val)  {drawFlag    = val;}; 
       
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
    std::vector<G4double> Evert;
    std::vector<G4int>    Nvert;
    
    Tst26EventMessenger* eventMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
