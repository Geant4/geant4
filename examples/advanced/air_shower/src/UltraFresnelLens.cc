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
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraFresnelLens.cc
//    ****************************************************
//
//    Class for definition of the Ultra Fresnel Lens. 
//    An  UltraFresnelLens object is created in the UltraDetectorConstruction class.
//    This class makes use of the UltraFresnelLensParameterisation class.
//    The lens profile is define through the GetSagita method.   
//
#include <cmath>
#include "UltraFresnelLens.hh"
#include "UltraFresnelLensParameterisation.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

UltraFresnelLens::UltraFresnelLens( G4double Diameter, G4int nGrooves, G4Material* Material, G4VPhysicalVolume * MotherPV, G4ThreeVector Pos)
{

 LensMaterial    = Material ;
 LensDiameter    = Diameter ;
 NumberOfGrooves = nGrooves;
 GrooveWidth     = (Diameter/2.0)/nGrooves ;
 LensPosition    = Pos ;

 if( GrooveWidth <= 0 ){
   G4Exception("UltraFresnelLens::UltraFresnelLens()","AirSh001",FatalException,
	       "UltraFresnelLens constructor: GrooveWidth<=0");
   }


G4double  Rmin1 = (NumberOfGrooves-1)*(GrooveWidth) ;
G4double  Rmax1 = (NumberOfGrooves-0)*(GrooveWidth) ;
LensThickness   = GetSagita(Rmax1)-GetSagita(Rmin1) ; // Height of the highest groove 

BuildLens(MotherPV) ;

}


UltraFresnelLens::~UltraFresnelLens(  )
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UltraFresnelLens::BuildLens(G4VPhysicalVolume *MotherPV){

G4double StartPhi = 0.0 ;
G4double DeltaPhi = twopi ;

G4double  LensMotherDz = LensThickness ;


G4Tubs *LensMotherCylinder 
    = new G4Tubs("LensMotherCylinder",0.0*mm,LensDiameter/2.0,LensMotherDz/2.0,StartPhi,DeltaPhi);

G4Material *LensMotherMaterial = MotherPV->GetLogicalVolume()->GetMaterial() ;
G4LogicalVolume *LensMotherLV 
    = new G4LogicalVolume(LensMotherCylinder,LensMotherMaterial,"LensMotherLV",0,0,0);
G4VPhysicalVolume *LensMotherPV 
    = new G4PVPlacement(0,LensPosition,"LensMotherPV",LensMotherLV,MotherPV,false,0);

LensMotherLV->SetVisAttributes (G4VisAttributes::GetInvisible());


G4Cons *solidGroove 
  = new G4Cons("Groove",40.0*mm,50.0*mm,40.0*mm,40.001*mm,1.0*mm,StartPhi,DeltaPhi);

G4LogicalVolume *logicalGroove 
    = new G4LogicalVolume(solidGroove,LensMaterial,"Groove_log",0,0,0);

G4VPVParameterisation *FresnelLensParam = new UltraFresnelLensParameterisation(this);

LensPhysicalVolume = 
 new G4PVParameterised("LensPV",logicalGroove,LensMotherPV,kZAxis,NumberOfGrooves,FresnelLensParam) ;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double UltraFresnelLens::GetSagita(G4double radius) 
{

  G4double Conic = -1.0;
  G4double Curvature = 0.00437636761488/mm ;
  G4double Aspher[8] = {      4.206739256e-05/(mm),
                                9.6440152e-10/(mm3),
                               -1.4884317e-15/(mm2*mm3),
			                  0.0/(mm*mm3*mm3),
			                  0.0/(mm3*mm3*mm3),
			                  0.0/(mm2*mm3*mm3*mm3),
			                  0.0/(mm*mm3*mm3*mm3*mm3),
			                  0.0/(mm*3*mm3*mm3*mm3*mm3)
                            };

  G4double TotAspher = 0.0*mm ;

  for(G4int k=1;k<9;k++){
    TotAspher += Aspher[k-1]*std::pow(radius,2*k) ;
  }

  G4double ArgSqrt = 1.0-(1.0+Conic)*std::pow(Curvature,2)*std::pow(radius,2) ; 

  if (ArgSqrt < 0.0){
    G4Exception("UltraFresnelLens::GetSagita()","AirSh002",
		FatalException,
		"UltraFresnelLensParameterisation::Sagita: Square Root of <0 !");
  }
  G4double Sagita_value = Curvature*std::pow(radius,2)/(1.0+std::sqrt(ArgSqrt)) + TotAspher;

  return Sagita_value ;

                             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


