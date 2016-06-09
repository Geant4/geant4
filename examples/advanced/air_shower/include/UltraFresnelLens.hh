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
//   **********************************************
//   *        UltraFresnelLens.hh
//   **********************************************
//
//    Class for definition of the Ultra Fresnel Lens. 
//    An  UltraFresnelLens object is created in the UltraDetectorConstruction class.
//    This class makes use of the UltraFresnelLensParameterisation class.
//    The lens profile is define through the GetSagita method.   
//
#ifndef UltraFresnelLens_H
#define UltraFresnelLens_H

#include "globals.hh"
#include "G4ThreeVector.hh"
 
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;


class UltraFresnelLens{

public:

UltraFresnelLens( G4double Diameter, G4int nGrooves, G4Material* Material, G4VPhysicalVolume * MotherPV, G4ThreeVector Pos);

~UltraFresnelLens();

  inline G4VPhysicalVolume*  GetPhysicalVolume() {return LensPhysicalVolume ;}
  inline G4Material*         GetMaterial()       {return LensMaterial ;}
  inline G4double            GetDiameter()       {return LensDiameter ;}
  inline G4double            GetThickness()      {return LensThickness ;}
  inline G4double            GetGrooveWidth()    {return GrooveWidth ;}
  inline G4int               GetNumberOfGrooves(){return NumberOfGrooves ;} 

         G4double            GetSagita(G4double) ;


private:

void     BuildLens(G4VPhysicalVolume *) ;

G4double      LensDiameter;
G4int         NumberOfGrooves;
G4Material*   LensMaterial;
G4ThreeVector LensPosition ;

G4double      GrooveWidth ;
G4double      LensThickness ;


G4VPhysicalVolume* LensPhysicalVolume ;

};

#endif 
