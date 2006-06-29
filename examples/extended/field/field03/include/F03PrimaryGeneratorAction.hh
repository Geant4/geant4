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
// $Id: F03PrimaryGeneratorAction.hh,v 1.4 2006-06-29 17:18:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F03PrimaryGeneratorAction_h
#define F03PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class F03DetectorConstruction;
class F03PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F03PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    F03PrimaryGeneratorAction(F03DetectorConstruction*);    
   ~F03PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);
    void SetRndmFlag(G4String val) { rndmFlag = val;}
    void Setxvertex(G4double x) ;
    void Setyvertex(G4double y) ;
    void Setzvertex(G4double z) ;

    static G4String GetPrimaryName() ;                

  private:
    G4ParticleGun*                particleGun;	// pointer a to G4 service class
    F03DetectorConstruction*      F03Detector;  // pointer to the geometry
    
    F03PrimaryGeneratorMessenger* gunMessenger; // messenger of this class
    G4String                      rndmFlag;     // flag for random impact point       

    static G4String thePrimaryParticleName ;
    G4double xvertex,yvertex,zvertex;
    G4bool vertexdefined ;

};

#endif


