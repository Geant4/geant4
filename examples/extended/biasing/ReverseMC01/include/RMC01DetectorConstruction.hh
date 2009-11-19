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
// $Id: RMC01DetectorConstruction.hh,v 1.1 2009-11-19 22:41:18 ldesorgh Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//////////////////////////////////////////////////////////////
//      Class Name:	RMC01DetectorConstruction
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
//////////////////////////////////////////////////////////////
// CHANGE HISTORY
//--------------
//      ChangeHistory:
//	 	17-11-2009 creation by L. Desorgher
//
//-------------------------------------------------------------
#ifndef RMC01DetectorConstruction_h
#define RMC01DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class RMC01DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RMC01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    RMC01DetectorConstruction();       
   ~RMC01DetectorConstruction();

  public:
    G4VPhysicalVolume* Construct();
    void UpdateGeometry();
    
    void SetSensitiveVolumeHeight(G4double h);
    void SetSensitiveVolumeRadius(G4double r);
    
    void SetShieldingThickness(G4double d);
    
  private:
     
    RMC01DetectorMessenger* detectorMessenger;  //pointer to the Messenger
    
    //Geometrical parameters
    //---------------------- 
     G4double shield_Thickness;
     G4double sensitive_cylinder_H;
     G4double sensitive_cylinder_Rout;
   
    
     
     
     
     
     
     
      
private:
    
     void DefineMaterials();
     
     G4VPhysicalVolume* ConstructSimpleGeometry(); 
     
};



#endif
