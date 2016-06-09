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
// ********************************************************************
// Code developed by:
//  S.Larsson
//
//    ********************************************
//    *                                          *
//    *    PurgMagPrimaryGeneratorAction.hh     *
//    *                                          *
//    ********************************************
//
// $Id: PurgMagPrimaryGeneratorAction.hh,v 1.2 2004/06/18 09:17:47 gunter Exp $
// GEANT4 tag $Name: geant4-06-02 $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef PurgMagPrimaryGeneratorAction_h
#define PurgMagPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class PurgMagDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class PurgMagPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PurgMagPrimaryGeneratorAction(PurgMagDetectorConstruction*);    
  ~PurgMagPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event*);
  void SetRndmVertex(G4bool val) { rndmVertex = val;} 
  
private:
  G4ParticleGun*                  particleGun;
  PurgMagDetectorConstruction*    PurgMagDetector;
  
  G4bool                       rndmVertex;      
};

#endif


