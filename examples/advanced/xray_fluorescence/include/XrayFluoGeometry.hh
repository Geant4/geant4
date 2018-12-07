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
//
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  22 Aug 2003  Alfonso Mantero Created 
//
// -------------------------------------------------------------------

#ifndef XrayFluoGeometry_hh
#define XrayFluoGeometry_hh 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4Material;
class XrayFluoVDetectorType;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoGeometry: public G4VUserDetectorConstruction
{

public:
  
  ~XrayFluoGeometry();
  
  virtual G4VPhysicalVolume* Construct()=0;
  
  virtual void UpdateGeometry()=0;


  virtual void SetSampleMaterial(G4String)=0;

  virtual void SetDetectorType(G4String)=0;

  //virtual static XrayFluoDetectorConstruction* GetInstance();

  virtual void PrintApparateParameters()=0; 

  virtual XrayFluoVDetectorType* GetDetectorType()=0;


  virtual G4double GetWorldSizeZ()=0;
  virtual G4double GetWorldSizeXY()=0;
  
  virtual G4Material* GetSampleMaterial()=0;

protected:
  
  XrayFluoGeometry();

  
// private:
  
//   virtual void DefineDefaultMaterials()=0;
//   virtual G4VPhysicalVolume* ConstructApparate()=0;

//   //calculates some quantities used to construct geometry
//   virtual void ComputeApparateParameters()=0;

};

#endif






