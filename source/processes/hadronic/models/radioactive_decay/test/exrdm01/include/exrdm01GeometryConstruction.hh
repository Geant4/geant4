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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef exrdm01GeometryConstruction1_h
#define exrdm01GeometryConstruction1_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		exrdm01GeometryConstruction.hh
//
// Version:		0.b.3
// Date:		29/02/00
// Author:		P R Truscott
// Organisation:	DERA UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		12115/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
// This class is used to define the geometry and materials to which the SSAT
// is to be applied.  Any GEANT4-compatible geometry description may be used.
// This example considers a non-concentric sphere of radius 30 cm displaced
// in the z-direction by 7 cm containing aluminium.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
// exrdm01GeometryConstruction ()
//    Constructor which defines the size of the experiment hall and the spheres.
//
// ~exrdm01GeometryConstruction ()
//    Destructor.
//
// G4VPhysicalVolume *construct ()
//    Defines the materials and logical/physical volumes making up the geometry.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29/02/00, P R Truscott, DERA UK
// 0.b.3 release.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
///////////////////////////////////////////////////////////////////////////////
//
class exrdm01GeometryConstruction : public G4VUserDetectorConstruction
{
  public:
    exrdm01GeometryConstruction ();
    ~exrdm01GeometryConstruction ();

  public:
     G4VPhysicalVolume *Construct ();

  private:
     G4int nel;

     G4double universe_x;
     G4double universe_y;
     G4double universe_z;

     G4double aSphere_r1;
     G4double aSphere_r2;
};
///////////////////////////////////////////////////////////////////////////////
#endif
