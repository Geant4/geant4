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
// $Id: XrayFluoGeometry.hh
// GEANT4 tag $Name: xray_fluo-V03-02-00
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






