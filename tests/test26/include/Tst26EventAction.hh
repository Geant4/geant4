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
// $Id: Tst26EventAction.hh,v 1.6 2006-06-29 21:53:09 gunter Exp $
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

    
