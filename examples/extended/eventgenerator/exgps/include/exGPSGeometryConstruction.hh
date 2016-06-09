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
#ifndef exGPSGeometryConstruction1_h
#define exGPSGeometryConstruction1_h 1
//
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;

////////////////////////////////////////////////////////////////////////////////
//
class exGPSGeometryConstruction : public G4VUserDetectorConstruction
{
  public:
    exGPSGeometryConstruction ();
    ~exGPSGeometryConstruction ();

  public:
     G4VPhysicalVolume *Construct ();

  private:

  G4VPhysicalVolume * universe_phys ;
  G4VPhysicalVolume * Al_phys ;
  G4VPhysicalVolume * aSphere_phys;

};
////////////////////////////////////////////////////////////////////////////////
#endif

