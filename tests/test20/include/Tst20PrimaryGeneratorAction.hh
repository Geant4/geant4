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
// $Id: Tst20PrimaryGeneratorAction.hh,v 1.3 2001-07-11 10:10:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst20PrimaryGeneratorAction_h
#define Tst20PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class Tst20DetectorConstruction;
//class Tst20PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst20PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Tst20PrimaryGeneratorAction(Tst20DetectorConstruction*);    
   ~Tst20PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);
    void SetRndmFlag(G4String val) { rndmFlag = val;}
    void Setxvertex(G4double x) ;
    void Setyvertex(G4double y) ;
    void Setzvertex(G4double z) ;

    static G4String GetPrimaryName() ;                

  private:
    G4ParticleGun*                particleGun;	//pointer a to G4 service class
    Tst20DetectorConstruction*      Tst20Detector; //pointer to the geometry
    
  //  Tst20PrimaryGeneratorMessenger* gunMessenger; //messenger of this class
    G4String                      rndmFlag;	//flag for a random impact point       

    static G4String thePrimaryParticleName ;
    G4double xvertex,yvertex,zvertex;
    G4bool vertexdefined ;

};

#endif


