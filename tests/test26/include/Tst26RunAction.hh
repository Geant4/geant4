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
// $Id: Tst26RunAction.hh,v 1.4 2006-06-29 21:53:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#ifndef Tst26RunAction_h
#define Tst26RunAction_h 1

#include "G4UserRunAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Tst26PrimaryGeneratorAction;
class G4Run;


class Tst26RunAction : public G4UserRunAction
{
  public:
  
    Tst26RunAction(Tst26PrimaryGeneratorAction*);
   ~Tst26RunAction();

    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);
    
    void AddParticle(G4int, G4int);
    void AddEvent(G4double, G4double, G4double, G4double, 
                  G4double, G4double, G4double, G4int);
 
  private:
  
    Tst26PrimaryGeneratorAction* Tst26Kin;   
    G4int evtNo;
    G4int calNo;
    G4double beamE;
    G4double trancFactor;
    G4double e[8];
    G4double e2[8];
    G4double s[8];
    G4int sgam;
    G4int sgamV;
    G4int sgamM;
    G4int sel;
    G4int selV;
    G4int selM;
    G4int spos;
    G4int sposV;
    G4int sposM;
};

#endif

