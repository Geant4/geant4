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
