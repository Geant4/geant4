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
// $Id: GammaRayTelPrimaryGeneratorAction.hh,v 1.4 2001-07-11 09:56:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelPrimaryGeneratorAction  ------
//           by G.Santin, F.Longo & R.Giannitrapani (30 nov 2000)
//
// ************************************************************



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelPrimaryGeneratorAction_h
#define GammaRayTelPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class GammaRayTelDetectorConstruction;
class GammaRayTelPrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  GammaRayTelPrimaryGeneratorAction(GammaRayTelDetectorConstruction*);    
  ~GammaRayTelPrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event*);
  void SetRndmFlag(G4String val) { rndmFlag = val;}
  void SetSourceType(G4int val) { nSourceType = val;}
  void SetSpectrumType(G4int val) { nSpectrumType = val;}
  void SetVertexRadius(G4double val) { dVertexRadius = val;}
  
private:
  G4ParticleGun*                particleGun;	  
  GammaRayTelDetectorConstruction*    GammaRayTelDetector;  
  GammaRayTelPrimaryGeneratorMessenger* gunMessenger; 
  G4String                      rndmFlag;    //flag for a random impact point
  G4int                         nSourceType;
  G4double                      dVertexRadius;
  G4int                         nSpectrumType;
};

#endif



