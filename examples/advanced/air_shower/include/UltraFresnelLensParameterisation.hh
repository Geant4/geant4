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
//   *        UltraFresnelLensParameterisation.hh
//   **********************************************
//
//    Class derived from G4VPVParameterisation and used to define a Fresnel lens geometry
//    through a parameterised replication of G4Cons volumes. These volumes are frustra 
//    of cones describing the lens grooves.
//    An  UltraFresnelLensParameterisation object is created in the UltraFresnelLens class
//
#ifndef UltraFresnelLensParameterisation_H
#define UltraFresnelLensParameterisation_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"

class UltraFresnelLens ;

class G4VPhysicalVolume;
class G4Cons;

// Dummy declarations to get rid of warnings ...
class G4Box;
class G4Trd;
class G4Trap;
class G4Sphere;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Orb;
class G4Polycone;
class G4Polyhedra;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class UltraFresnelLensParameterisation : public G4VPVParameterisation
{ 
  public:
  
    UltraFresnelLensParameterisation(UltraFresnelLens* ) ;

    virtual				 
   ~UltraFresnelLensParameterisation();
   
    void ComputeTransformation (const G4int GrooveNo,
                                G4VPhysicalVolume* physVol) const;
    
    void ComputeDimensions (G4Cons & Groove, const G4int GrooveNo,
                            const G4VPhysicalVolume* physVol) const;
private:

  private:  // Dummy declarations to get rid of warnings ...
    void ComputeDimensions (G4Box&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Tubs&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trd&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,const G4VPhysicalVolume*) const {}

  private:

    UltraFresnelLens* FresnelLens ;
    G4int NumberOfGrooves ;
    G4double GrooveWidth ;
    G4double dZOffset ;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


