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
// $Id: FCALPrimaryGeneratorAction.hh,v 1.5 2004/11/29 18:03:05 ribon Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FCALPrimaryGeneratorAction_h
#define FCALPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class FCALTestbeamSetup;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FCALPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    FCALPrimaryGeneratorAction();    
   ~FCALPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);
  //    void SetRndmFlag(G4String val) { rndmFlag = val;}

  private:
    G4ParticleGun*                particleGun;	  //pointer a to G4 service class
    
  
  private:
  G4int Ievent;
  G4double* X;
  G4double* Y;
  G4double* Z;
  G4double* Cos_X;
  G4double* Cos_Y;
  G4double* Cos_Z;

  G4int Nevent;

};

#endif


