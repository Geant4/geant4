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
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraFresnelLensParameterisation.cc
//    ****************************************************
//
//    Class derived from G4VPVParameterisation and used to define a Fresnel lens geometry
//    through a parameterised replication of G4Cons volumes. These volumes are frustra 
//    of cones describing the lens grooves.
//    An  UltraFresnelLensParameterisation object is created in the UltraFresnelLens class
//
#include <cmath>
#include "UltraFresnelLensParameterisation.hh"
#include "UltraFresnelLens.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Cons.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UltraFresnelLensParameterisation::UltraFresnelLensParameterisation(UltraFresnelLens* Lens)  
{
   
   FresnelLens     = Lens ;
   GrooveWidth     = Lens->GetGrooveWidth() ;
   NumberOfGrooves = Lens->GetNumberOfGrooves() ;

   dZOffset = Lens->GetSagita((NumberOfGrooves-0)*(GrooveWidth)) - 
              Lens->GetSagita((NumberOfGrooves-1)*(GrooveWidth)) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UltraFresnelLensParameterisation::~UltraFresnelLensParameterisation()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UltraFresnelLensParameterisation::ComputeTransformation
(const G4int GrooveNo, G4VPhysicalVolume* physVol) const
{

  G4double  Rmin1 = (GrooveNo+0)*(GrooveWidth) ;
  G4double  Rmax1 = (GrooveNo+1)*(GrooveWidth) ;

  G4double     dZ = FresnelLens->GetSagita(Rmax1) - FresnelLens->GetSagita(Rmin1) ;

  if (dZ <= 0.0){
     G4Exception("UltraFresnelLensParameterisation::ComputeTransformation: Groove depth<0 !");
  }
 
  G4ThreeVector origin(0,0,(dZ-dZOffset)/2.);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UltraFresnelLensParameterisation::ComputeDimensions
(G4Cons& Groove, const G4int GrooveNo, const G4VPhysicalVolume*) const         
{
  G4double  Rmin1 = (GrooveNo+0)*(GrooveWidth) ;
  G4double  Rmax1 = (GrooveNo+1)*(GrooveWidth) ;

  G4double  Rmin2 = Rmin1 ;
  G4double  Rmax2 = Rmin2+0.0001*mm ;
  
  G4double  dZ = FresnelLens->GetSagita(Rmax1) - FresnelLens->GetSagita(Rmin1) ;

  if (dZ <= 0.0){
     G4Exception("UltraFresnelLensParameterisation::ComputeDimensions: Groove depth<0 !");
  }


  Groove.SetInnerRadiusMinusZ(Rmin1) ;
  Groove.SetOuterRadiusMinusZ(Rmax1) ;

  Groove.SetInnerRadiusPlusZ(Rmin2) ;
  Groove.SetOuterRadiusPlusZ(Rmax2) ;

  Groove.SetZHalfLength(dZ/2.) ;

#ifdef ULTRA_VERBOSE

G4cout << "UltraFresnelLensParameterisation: GrooveNo " << GrooveNo+1 << 
" Rmin1, Rmax1(mm): " << Rmin1/mm <<" "<<  Rmax1/mm << " dZ(mm) " << dZ/mm << G4endl ;
#endif



}

