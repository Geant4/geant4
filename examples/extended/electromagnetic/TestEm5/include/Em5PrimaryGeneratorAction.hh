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
// $Id: Em5PrimaryGeneratorAction.hh,v 1.2.4.1 2001/06/28 19:07:05 gunter Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em5PrimaryGeneratorAction_h
#define Em5PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class Em5DetectorConstruction;
class Em5PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em5PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Em5PrimaryGeneratorAction(Em5DetectorConstruction*);    
   ~Em5PrimaryGeneratorAction();

  public:
    void SetDefaultKinematic();   
    void GeneratePrimaries(G4Event*);
    static G4String GetPrimaryName() ;                

  private:
    G4ParticleGun*                particleGun;
    Em5DetectorConstruction*      Em5Detector;
    
    static G4String thePrimaryParticleName;
    
    Em5PrimaryGeneratorMessenger* gunMessenger;     
};

#endif


