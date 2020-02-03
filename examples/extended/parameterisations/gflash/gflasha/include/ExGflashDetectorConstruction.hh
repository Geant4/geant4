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
/// \file ExGflashDetectorConstruction.hh
/// \brief Definition of the ExGflashDetectorConstruction class
//
#ifndef ExGflashDetectorConstruction_h
#define ExGflashDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "ExGflashSensitiveDetector.hh"
#include "G4ThreeVector.hh"
#include "G4Cache.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4Region;

class GFlashHomoShowerParameterisation;
class GFlashShowerModel;
class GFlashHitMaker;
class GFlashParticleBounds;
class ExGflashMessenger;

class ExGflashDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  ExGflashDetectorConstruction();
  ~ExGflashDetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  
  void SetLBining (G4ThreeVector);
  void SetRBining (G4ThreeVector);
  void SetVerbose(G4int val)  {fVerbose = val;}

  void SetMaterial(G4String mat);

  G4int    GetVerbose() const    {return fVerbose;}

  G4int    GetnLtot() const      {return fNLtot;}
  G4int    GetnRtot() const      {return fNRtot;}
  G4double GetdLradl() const     {return fDLradl;}
  G4double GetdRradl() const     {return fDRradl;}

  G4double GetSDRadLen() const   {return fSDRadLen;}
  
private:
  G4int fNbOfCrystals; // cube of nb x nb crystals

  G4double fCrystalWidth; // x,y size
  G4double fCrystalLength;// z size

  G4LogicalVolume*    fCrystal_log;
  G4Material* fDetMat;
  G4Region*           fRegion;

  G4double fSDRadLen;  // SD material Rad Lenght

  G4int fVerbose;
  G4int    fNLtot,    fNRtot;       // nb of bins: longitudinal and radial
  G4double fDLradl,   fDRradl;      // bin thickness (in radl unit)

  ExGflashMessenger*  fGflashMessenger;

  static G4ThreadLocal GFlashShowerModel* fFastShowerModel; 
  static G4ThreadLocal GFlashHomoShowerParameterisation* fParameterisation;
  static G4ThreadLocal GFlashParticleBounds* fParticleBounds;
  static G4ThreadLocal GFlashHitMaker* fHitMaker;
};

#endif
