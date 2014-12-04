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

#include "TstTarget.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
//#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4RunManager.hh"

// #include <fstream>
#include <string>
#include <iostream>
#include <sstream>
// #include <iomanip>

TstTarget::TstTarget()
   : G4VUserDetectorConstruction(),
     fX1(100.), fX2(100.), fX3(100.), // set default size
     fShape("G4Box"),
     fMatName(""),
     fMaterial(0),
     fWorld(0), fSubWorld(0),
     fLogTarget(0), fPhysTarget(0)
{
}

TstTarget::~TstTarget()
{

   // Do I need to delete log/phys volumes ?
   // Or are they taken care of by someone else ?

}

G4Material* TstTarget::ResetMaterial( G4String mat )
{
   

   if ( mat.find("G4") != std::string::npos ) 
   {
      fMatName = mat ;
      fMaterial = G4NistManager::Instance()->FindOrBuildMaterial(fMatName);
      if (!fMaterial) 
      {
         G4cout << "Material <" << fMatName << "> is not found" << G4endl;
         exit(1);
       }
      return fMaterial;
   }

   std::ostringstream osMat(std::ios_base::out|std::ios_base::app);
   osMat.clear();
   osMat.str("G4_");
   osMat << mat;
   fMatName = osMat.str();

   fMaterial = G4NistManager::Instance()->FindOrBuildMaterial(fMatName);
   if (!fMaterial) 
   {
      G4cout << "Material <" << fMatName << "> is not found" << G4endl;
      exit(1);
   }
      
   return fMaterial;

}

void TstTarget::ResetMaterial( G4Material* mat )
{

   if ( !mat )
   {
      G4cout << "Invalid Material" << G4endl;
      exit(1);
   }
   
   
   if ( !mat->GetName() )
   {
      G4cout << "Invalid Material" << G4endl;
      exit(1);      
   }
    
   fMaterial = mat;
   fMatName = mat->GetName();

   return;

}

void TstTarget::ResetGeom()
{

   // Clean old geometry, if any.
   //
   G4GeometryManager::GetInstance()->OpenGeometry();
   G4PhysicalVolumeStore::GetInstance()->Clean();
   G4LogicalVolumeStore::GetInstance()->Clean();
   G4SolidStore::GetInstance()->Clean();

   return;

}

G4VPhysicalVolume* TstTarget::Construct()
{

   ResetGeom();
   
   G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
   
   // G4Box* solidWorld = new G4Box( "world_solid", 10.*m, 10.*m, 10.*m );
   fWorld = new G4PVPlacement( 0, G4ThreeVector(), "world_phys",
                               new G4LogicalVolume( new G4Box( "world_solid", 10.*m, 10.*m, 10.*m ),
			                            mat, "world_log", 0, 0, 0 ), // mat=Vacuum !!!
			       0,
			       false, 0 );
   
   // G4Box* solidSubWorld = new G4Box( "subworld_solid", 9.999*m, 9.999*m, 9.999*m );
   fSubWorld = new G4PVPlacement( 0, G4ThreeVector(), "subworld_phys",
                                  new G4LogicalVolume( new G4Box( "subworld_solid", 9.999*m, 9.999*m, 9.999*m ),
				                       mat, "subworld_log", 0, 0, 0 ), // mat=Vacuum again
				  fWorld,
				  false, 0 );

   ConstructTarget();
   
   return fWorld;

}

G4VPhysicalVolume* TstTarget::ConstructTarget()
{

  // Clean old geometry, if any.
  //
/*
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
*/

  G4CSGSolid* solid = 0;
  
  if ( fShape == "G4Box" )
  {
     solid = new G4Box ( "target_solid", 0.5*fX1, 0.5*fX2, 0.5*fX3 ); // size is already done in the right units (in TstReader)
  }
  else if ( fShape = "G4Tubs" )
  {
     solid = new G4Tubs( "solid", fX1, fX2, 0.5*fX3, 0., 2.*pi );
  }
  
  fLogTarget = new G4LogicalVolume( solid, fMaterial, "target_log", 0, 0, 0);
  fPhysTarget = new G4PVPlacement( 0, G4ThreeVector(), "target_phys",
                                   fLogTarget,
				   fSubWorld, // mother physical volume, which can be NULL
				   false, 0 );

  return fPhysTarget; 

}
