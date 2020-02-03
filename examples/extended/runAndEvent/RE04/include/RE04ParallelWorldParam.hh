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
/// \file runAndEvent/RE04/include/RE04ParallelWorldParam.hh
/// \brief Definition of the RE04ParallelWorldParam class
//
//
#ifndef RE04ParallelWorldParam_h
#define RE04ParallelWorldParam_h 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"

class G4VPhysicalVolume;
class G4Material;

// Dummy declarations to get rid of warnings ...
class G4Box;
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Ellipsoid;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;

//
/// Parameterisation class for volumes in a parallel world 
///
/// - void ComputeTransformation(const G4int copyNo,
///                              G4VPhysicalVolume *physVol) const
///     returns a transformation with the physical volume of the 2nd argument
///     according to copyNo.
///       copyNo = 0 : the volume is placed at (-10*cm,-10*cm,0*cm)
///              = 1 : the volume is placed at (10*cm,10*cm,0*cm)
/// - G4Material* ComputeMaterial(const G4int copyNo,
///                               G4VPhysicalVolume* currentVol,
///                               const G4VTouchable* parentTouch=0);
///     returns a material according to copyNo.
///       copyNo = 0 : water
///              = 1 : lead
//
class RE04ParallelWorldParam : public G4VPVParameterisation
{
  public:

    RE04ParallelWorldParam();
   ~RE04ParallelWorldParam();

  public:

    void ComputeTransformation(const G4int copyNo,
                                     G4VPhysicalVolume *physVol) const;
    G4Material* ComputeMaterial(const G4int copyNo,
                                      G4VPhysicalVolume* currentVol,
                                const G4VTouchable* parentTouch=0);

  public:
    virtual void ComputeDimensions (G4Box&,const G4int,const G4VPhysicalVolume*)
      const {}

  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const {}

  private:
    G4Material* fWater;
    G4Material* fPb;
};

#endif

